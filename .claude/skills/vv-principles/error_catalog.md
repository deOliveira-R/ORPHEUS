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
`docs/theory/methods/sn/index.rst` §Investigation History.

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
  `docs/theory/methods/diffusion_1d.rst` "Investigation history".

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

**L0 test that catches it (post-closure evidence):**
`tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_sphere`
— parametrized over ``inner_solver ∈ {"source_iteration", "krylov"}``,
``n_cells ∈ {20, 40, 80}``, ``n_ord ∈ {8, 16}``.  Both inner solvers must
reach the closed-form analytical streaming-equilibrium answer
``φ = Q/(Σ_t(1-c))`` with per-ordinate ``ψ_n = φ/Σw`` at ``rtol = 1e-9``.
Carries ``@pytest.mark.catches("ERR-026", "ERR-048")`` since the
PR-CLEANUP-CODE retirement of the historical evidence ledger
``tests/sn/test_sweep_operator_inconsistency.py`` (2026-05-16).  The
canonical SI-vs-Krylov three-way standoff that documented the WDD
fixed-point bias has been superseded by this L0 test (which now PASSES
on both solvers — the cleanest possible evidence of ERR-026 closure).

Additional ERR-026-tagged tests:
``tests/sn/test_phase_c_gates.py``, ``tests/sn/test_phase_c_mms.py``,
``tests/sn/spatial/test_psi_half_angle_seed.py``,
``tests/sn/spatial/test_sweep_vs_apply_consistency.py``,
``tests/sn/spatial/test_apply_matvec_cylinder_invariants.py`` —
collectively the post-closure regression net.

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
:doc:`docs/theory/methods/sn/index`.

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

What Wave H Phase C added (commits `eae6f05`..., GH #168 Phase C):

Phase C (2026-05-12) shipped the **sweep-frame apply matvec rewrite**
that aligns the apply matvec's spatial closure with the sweep's WDD
form. The architectural changes:

1. `transport_operator_matvec_spherical` and `_cylindrical`
   rewritten as one sweep iteration semantically — WDD diamond
   `psi_face_out = 2*psi_cell - psi_face_in` cell-by-cell along
   the DAG-ordered cell sequence
   (`SNMesh.iter_cells_by_direction(±1)`, added in Phase C); BC
   trace law applied ONCE at the boundary edge on the
   WDD-propagated outflow face values (the §16A.3 contract).
   Vectorised across ordinates via outgoing_mask / incoming_mask.

2. `solution_to_angular_flux_spherical` (+ cylindrical alias)
   simplified to return only `fi` (the Phase A
   `(fi, boundary_face_flux)` companion array retires).

3. Phase A's `BoundaryFaceFlux` Protocol RETIRED entirely:
   `orpheus/sn/spatial/boundary_face_flux.py` deleted, the 21
   foundation tests at `tests/sn/spatial/test_boundary_face_flux.py`
   deleted, `SNMesh.boundary_face_flux` field removed,
   `boundary_face_flux_closure` kwarg removed from the matvec
   signatures.

Empirical Gate 1.1 outcome (the load-bearing decision point):

Phase C's Gate 1.1 probe (per-ordinate flat-flux residual on
reflective curvilinear) with `MorelMontryAngularSweep` as the
candidate default FAILED on sphere. Cylindrical-MMS passes
(per-level α-telescoping happens to align with the WDD pole-face
initial condition `ψ_face_in(pole) = ψ_cell[0]` used in the
sweep-frame rewrite); spherical-MMS does NOT (the pole-face
initial condition interacts with the canonical Hébert §3.9.4
angular recurrence in a way that breaks the per-ordinate flat-flux
invariant). Per user constraint 7, the curvilinear default stays
`LegacyTauSymmetricInterpolation`; the default flip is **DEFERRED
to Phase D**.

What Phase C ships (architectural):

- Sweep-frame matvec rewrite with WDD spatial closure
  cell-by-cell across the DAG-ordered cell sequence.
- BC trace law applied at the boundary edge on the WDD-propagated
  outflow face — honours the §16A.3 contract structurally.
- Phase A BoundaryFaceFlux Protocol retired entirely.
- All 6 curvilinear regression snapshots regenerated against the
  Phase C operator (commit `7497cec` completed the
  `sphere_2g_p1_aniso` regen after the catalog narrative landed in
  `d445a8f`).

What Phase C leaves open (Phase D scope):

The four `xfail-strict` curvilinear MMS tripwires STAY xfail:

- `tests/sn/test_mms_curvilinear.py::test_sn_spherical_mms_converges_second_order`
- `tests/sn/test_mms_curvilinear.py::test_sn_cylindrical_mms_converges_second_order`
- `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py::test_sn_spherical_aniso_mms_converges_second_order`
- `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py::test_sn_cylindrical_aniso_mms_converges_second_order`

These will close in Phase D once the pole-face spatial closure
refinement lets MorelMontryAngularSweep become the natural default.

Tests added in Phase C (`tests/sn/test_phase_c_*.py`):

- `test_phase_c_gates.py` — Gate 1.1 (parametrised over 3 pole
  closures × 2 Σ_t × 2 geometries; sphere-MMS marked xfail per
  the empirical finding), Gate 1.2 (apply determinism via
  `np.array_equal`), Gate 1.3 (apply ↔ apply_transpose
  reciprocity), Gate 1.4 (linearity), Gate 1.5 (BC trace
  contract).
- `test_phase_c_mms.py` — Gate 3.1 / 3.2 (spatial MMS, both
  marked xfail with ERR-026 PARTIAL CLOSURE rationale) + Gate 3.3
  (angular convergence at fixed nx — passes).
- `test_phase_c_crosscheck.py` — Gate 4.1 (k_∞ recovery PASSES)
  and Gate 4.2 placeholder (SKIP for Phase D follow-up).

Phase C ships ERR-026 status: still **PARTIAL CLOSURE** (the
sweep-frame matvec aligns the spatial-closure architecture — the
Phase B load-bearing precondition — but the angular default flip
+ pole-face spatial closure refinement is Phase D's scope).

What Wave H Phase D added (commits `9512459`..`c44fe9b`, GH #168 +
#192, 2026-05-12):

Phase D (2026-05-12) shipped the canonical Hébert §3.9.4 Eqs.
(3.432)-(3.435) **Carlson coupled-pole inward μ=−1 sweep** as the
half-angle face flux seed for the M-M angular recurrence. Architecture:

1. New module `orpheus/sn/spatial/psi_half_angle_seed.py` introduces
   the `PsiHalfAngleSeed` Protocol family with two strategies:
   - `ZeroSeed` (key=`"zero"`): reproduces Phase B's hardcoded
     `psi_half_left = 0` behaviour bit-identically (regression-safety
     ablation knob).
   - `CarlsonInwardSweep` (key=`"carlson_inward_sweep"`): canonical
     Hébert §3.9.4 inward zero-weight sweep at μ=−1. Returns the
     per-cell half-angle profile `(ng, nx)` the M-M α-cascade
     consumes via the redistribution coefficient.

2. `MorelMontryAngularSweep` gains a `psi_half_seed: PsiHalfAngleSeed
   = field(default_factory=CarlsonInwardSweep)` dataclass field
   (Option α — composition into the M-M closure rather than a
   sibling Protocol on SNMesh; the seed is M-M-specific because
   Legacy and Bailey closures have no `psi_half_left` variable).

3. `transport_operator_matvec_spherical` and `_cylindrical` build a
   `CarlsonSweepContext` bundling `(sigma_t, dr, mu_quad, weights,
   bc_outer_value)` and pass it to `pole_angular_closure` via an
   optional `carlson_context` kwarg. Legacy and BFF ignore the
   kwarg; M-M consumes it via its `psi_half_seed` strategy.

4. **Default flips**: `SNMesh.pole_angular_closure` default flipped
   `LegacyTauSymmetricInterpolation` → `MorelMontryAngularSweep`
   (which itself defaults to `CarlsonInwardSweep` for its seed).
   `solve_sn_fixed_source` curvilinear default `inner_solver` flipped
   `"source_iteration"` → `"krylov"`.

5. Gate 1.5 strengthened: BC trace contract test upgraded from
   `apply(0)=0` probe to capture-and-compare, asserting BOTH BC
   apply calls per matvec (Phase D Carlson context + Phase C trace
   law) match independently-reconstructed references at
   `rtol=atol=1e-14`. Parametrised over vacuum + reflective BCs.

The corrected injection-point (the load-bearing structural
finding):

Phase D's Step 2 diagnostic memo at
`.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md`
empirically established that the Phase D plan's original Protocol
injection point (`operator.py:734-740`'s `psi_face_in` for the WDD
outward sweep) is a **no-op on flat ψ** — the spatial pole-face IC
already coincides with the cell-centred value on flat ψ. The
load-bearing site is the M-M angular recurrence at
`orpheus/sn/spatial/pole_angular_closure.py:411`'s hardcoded
`psi_half_left = np.zeros(...)`. Replacing the zero seed with the
Carlson sweep's `φ̄_{1/2,i}` profile closes the per-ordinate
residual to ≤ 1e-15 on the Gate 1.1 sphere MMS probe.

Empirical evidence of partial closure:

- **Gate 1.1 sphere MMS (per-ordinate identity)**: 4 parametrised
  cases (sphere × Σ_t∈{0, 0.5} + cylinder × Σ_t∈{0, 0.5}) now
  XPASS under `MorelMontryAngularSweep` — residual collapses from
  O(10) (Phase C baseline, max |r| ≈ 18.88) to O(1e-15) (post-
  Phase-D). Identity-manifestation of ERR-026 CLOSED.
- **Convergence rate (L1 MMS magnitude)**: empirical probe on the
  L1 sphere isotropic ansatz at nx ∈ {20, 40, 80} under
  Krylov+Carlson gives L2 errors [6.39, 0.63, 0.11] with orders
  [3.33, 2.46] — PASSES O(h²). Under SI+Carlson the errors
  plateau at [0.083, 0.095, 0.098] with orders [-0.19, -0.04] —
  PLATEAUS (Krylov-flip is the correct default). Rate-manifestation
  of ERR-026 CLOSED.
- **Convergence magnitude (L1 MMS absolute bound)**: the four
  L1 MMS xfail-strict markers still XFAIL because the test's
  absolute-magnitude check `errors[-1] < 1e-3` fails at nx≤160
  (the Krylov+Carlson path's error at nx=80 is ~0.11; extrapolating
  with order 2.46 to nx=160 gives ~0.02 — still > 1e-3).
  Magnitude-manifestation of ERR-026 stays open: tracked as
  Issue #195.

Tests added or updated in Phase D (`tests/sn/spatial/test_psi_half_angle_seed.py`,
`tests/sn/test_phase_c_gates.py`, `tests/sn/test_snstreamingoperator.py`):

- `test_psi_half_angle_seed.py` (NEW) — 25 foundation + L0 + L1
  tests covering Protocol conformance, registry, immutability,
  shape contract, bit-identity for ZeroSeed, flat-ψ algebraic
  identity (reflective + varying C + vacuum nx=3 hand-calc),
  multi-region σ_t step, linearity for both seeds, structural
  independence (Carlson vs ZeroSeed on vacuum-BC probe), M-M
  default seed pinning.
- `test_phase_c_gates.py` — Gate 1.5 strengthened to capture-and-
  compare with positional discrimination of the two BC apply calls.
- `test_snstreamingoperator.py` — 3 tests updated: 1 docstring +
  threshold rewrite (now pins Phase D fix), 2 bit-identity tests
  threaded with `sn_mesh.pole_angular_closure`, 1 linearity
  tolerance relaxed `rtol=1e-13 → 1e-12` (principled per
  vv-principles §"Bit-identity vs principled-equivalence": new
  named intermediates, structurally-independent reference, drift
  bounded by ~25 × ULP).

Phase D status after Step 3 ships: ERR-026 **PARTIAL CLOSURE**
(narrowed scope — the per-ordinate identity AND convergence rate
manifestations are CLOSED; the L1 magnitude manifestation is
deferred to Issue #195 which will investigate whether the pre-
asymptotic transient is benign or hides a coefficient bug). Marker
removal blocked on #195; the marker reason strings updated to
attribute the failure to the magnitude check (not the rate or
identity, which Phase D restored).

Phase D follow-up issues:

- **#193** — BC-realizer level-locality invariant test (forward-
  looking risk for future tilted-BC kinds; cylindrical Carlson
  context commutativity assumption).
- **#194** — Wire `verifies('hebert-3-43X')` decorators on the L0
  algebraic identity tests OR remove orphan label declarations.
- **#195** — L1 curvilinear MMS magnitude check pre-asymptotic
  transient investigation + marker removal trigger.

What Wave H Phase E added (commits ``2d3e7f2``..``6708a4a``,
GH #168 Phase E, 2026-05-12):

Phase E focused on tightening the Gate 4.2 SN-vs-Variant-α cross-
check on heterogeneous MR by fixing a documented prototype TODO in
``orpheus/derivations/continuous/trajectory_resolvent/greens_function.py``.
The Variant α MR solvers (sphere + cylinder) had been deployed with
single-domain Gauss-Legendre quadrature on (0, R) — an inline
``follow-on improvement`` comment flagged it. Single-domain GL on a
piecewise-constant σ_t domain cannot resolve material kinks; this
manifested as non-monotone k_eff convergence under refinement
(n_r=24 → 7% off SN; n_r=32 → 22.7% off; n_r=64 → 14.7% off).

Phase E shipped composite per-region GL (helper
``_composite_per_region_gl``), wired into all three MR entry points.
Each region's GL is exact for polynomials of degree 2·n_k − 1 on
that region. Empirical post-Phase-E heterogeneous MR k_eff gap:
**1.99e-4 (sphere) / 1.75e-2 (cylinder)** at production quadrature,
down from 7–9%. Gate 4.2 rtol tightened 1e-1 → 2e-2 (sphere) /
3e-2 (cylinder). The pre-Phase-E ceilings on
``test_mr_interface_continuity_3region`` (cylinder) and
``test_garcia_case1_phi_matches_at_point[r=7.0]`` (sphere outer
surface) were re-calibrated since composite-GL nodes sit strictly
INSIDE regions (the spline must now extrapolate to the boundary,
trading interior accuracy for boundary artifact).

Phase E ALSO built a flux-shape comparison sentinel
(``test_phase_e_trajectory_resolvent_flux_shape_crosscheck``) which
discovered a separate finding: SN and Variant α agree on the
eigenVALUE (2e-4) but DISAGREE on the eigenVECTOR shape (65% per-
cell on sphere). Xfailed-strict pending Phase F investigation.

Phase E status: ERR-026 manifestation #4 (trajectory_resolvent
MR quadrature instability) **CLOSED**. NEW manifestation surfaced:
#6 (heterogeneous eigenvector shape) **OPEN** pending Phase F.

What Wave H Phase F added (commits ``<TBD>``,
GH #168 Phase F, 2026-05-12):

Phase F closed the structural twin of the Phase D bug: ``CarlsonInwardSweep``
(Hébert §3.9.4 Eqs. 3.432-3.435) was wired into the apply-matvec
path by Phase D Step 3, but the SI/sweep path
(``orpheus/sn/sweep.py::_sweep_1d_spherical`` and
``_sweep_1d_cylindrical``) still initialized the M-M angular
recurrence half-angle face flux to **zero**. That meant the snapshot
generator (which uses SI by default) produced fixed points under
the wrong seed, while the apply path (used by Krylov) produced
the right seed. The two paths solved DIFFERENT equations on the
same problem (n=40 sphere_2g_3reg: SI k_eff 1.38070, Krylov 1.38464,
2% disagreement diverging under refinement on per-cell shape).

Phase F factored a NEW helper
``orpheus/sn/spatial/psi_half_angle_seed.carlson_inward_sweep_from_source``
that runs the same Hébert inward sweep but driven by source
(``Q̄_i = Q_1d[i]/Σw``) instead of by the apply-path's current
``φ_0``. The SI/sweep dispatch sites at ``sweep.py:474`` (spherical,
``psi_angle = np.zeros((nx, ng))``) and ``sweep.py:634`` (cylindrical
per-level zero init) replaced with calls to this helper. The
apply-path math is unchanged (refactor at the helper level).

Empirical post-Phase-F:
- ``sf[0]/sf[1]`` on ``sphere_2g_3reg`` n=40: 0.522 → **0.778**;
  STABLE under refinement (was DIVERGING to 0.473 at n=320).
- ``cv(ψ@i=0)``: 0.520 → **0.404**; max/min(ψ) 6.4× → 1.16× (per-
  ordinate quasi-isotropy substantially restored, Pomraning 1989
  prediction approached).
- ``sf[-1]/sf[-2]`` (outer-reflective): 0.887 → **0.997** (outer-
  cell defect essentially CLOSED).
- SI-vs-Krylov k_eff gap converges **O(h)** (0.286% at n=40 →
  0.065% at n=160); both methods now converge to the same value.

NEW tests added:
- ``tests/sn/test_phase_c_gates.py::test_sweep_curvilinear_per_ordinate_flat_flux_residual``
  — Gate 1.6, the DUAL of Phase D Gate 1.1 for the SI/sweep path.
- ``tests/sn/spatial/test_sweep_vs_apply_consistency.py`` — 57
  foundation tests pinning apply-vs-sweep Carlson seed equivalence.

Phase F status: ERR-026 manifestation #6 (heterogeneous eigenvector
shape, structural pole-cell defect) **CLOSED**. The Phase E flux-
shape sentinel still xfails due to residual O(h) WDD spatial-closure
asymmetry between SI and Krylov paths; **reclassified as a new
manifestation #7 (convergence-rate of SI-vs-Krylov per-cell agreement)
PARTIAL** — tracked separately, the structural divergence is gone.

Updated ERR-026 manifestation table:

| # | Manifestation | Status |
|---|---|---|
| 1 | Spatial closure inconsistency | CLOSED by Phase C |
| 2 | Per-ordinate identity (apply-path Gate 1.1) | CLOSED by Phase D |
| 3 | Convergence rate (L1 MMS rate) | CLOSED by Phase D Krylov flip |
| 4 | trajectory_resolvent MR quadrature instability | CLOSED by Phase E |
| 5 | L1 absolute magnitude (errors[-1] < 1e-3) | **CLOSED by ERR-058 (#195)** — curvilinear isotropic MMS collapses O(h²) into the band (sphere 5.74e-5 @ nx=320, cyl 3.37e-5 @ nx=160) |
| 6 | Heterogeneous eigenvector shape (sweep-path Carlson seed) | **CLOSED by Phase F** (full closure via Phase G Step 2 ERR-048 pole-face IC + Carlson seed source fixes) |
| 7 | SI-vs-Krylov per-cell agreement (residual O(h) WDD asymmetry) | **CLOSED by ERR-058 (#195), verified + pinned by #196 (2026-06-12).** Phase G Step 2 (ERR-048) closed only the *L0 flat-field* twin-agreement (SI sweep ≡ apply-matvec on the streaming-equilibrium gauntlet); the *L1 heterogeneous eigenvalue* O(h) asymmetry this row names PERSISTED — hence #196 stayed open — because the shared closure seeds were still exact-on-flat / O(1)-wrong-on-non-flat (the ERR-058 defect). ERR-058's correct shared seeds (coupled-pole spatial + `AngularEdgeExtrapolation`) closed it at the eigenvalue level: SI ≡ Krylov to ~1.9e-11 k_eff / ~2.4e-10 flux-shape on sphere_2g_3reg n=40 (bug-era 3.9e-3 / ~30 %), pinned by `tests/sn/eigenvalue/test_keff_curvilinear.py::test_si_krylov_eigenvalue_equivalence_{sphere,cylinder}` |

Phase F follow-up (DONE 2026-06-12): the new issue tracking
manifestation #7 was filed as **#196** and is now **CLOSED**.
Neither Phase-F closure option was needed — option (a) "make SI
bit-identical to Krylov" and option (b) "flip the curvilinear
default to ``krylov``" both presupposed the sweep fixed point was
CORRECT and only the SI-vs-matvec arithmetic differed.  The
terminal diagnosis (ERR-058, #195) showed the shared fixed point
itself was WRONG on non-flat fields; fixing the closure seeds made
SI ≡ Krylov to the iteration floor on the SAME (now correct)
operator, so the curvilinear default stays ``source_iteration``
(reverted from the Phase-D Krylov flip — SI is ~10²× faster and now
bit-equivalent in fixed-source, floor-equivalent in eigenvalue).
The canary ``test_phase_e_trajectory_resolvent_flux_shape_crosscheck``
runs as a plain L1 test (xfail removed); the permanent
manifestation-#7 catcher is
``tests/sn/eigenvalue/test_keff_curvilinear.py::test_si_krylov_eigenvalue_equivalence_{sphere,cylinder}``.

**Anti-pattern surfaced by Phase F (proposed addition to
``vv-principles/SKILL.md``)**:

> Whenever a fix is applied to one of two structurally-mirrored
> production paths (apply-matvec vs SI/sweep, etc.), MUST audit
> the OTHER path for the same defect. Mode 3 wrong-term-
> initialization defects often appear in pairs; fixing one path
> without auditing its sister is a Cardinal Rule 2 (architecture)
> violation that ERR-026 instantiated twice.

### CLOSED 2026-06-12 (ERR-058 / Issue #195)

**Status: CLOSED.** The wrong-fixed-point family this entry opened
(2026-04-17, the 35–50% r=0 deviation that grows under refinement) is
terminally closed by the ERR-058 closure-seed fix — see the ERR-058
entry for the full mechanism: (a) the pole-face CELL-CENTRE spatial
seed → Carlson coupled-pole seed ψ(0,+μ) = ψ(0,−μ); (b) the M-M
half-angle thread's proxy-source Carlson seed (Q̄ = Σ_t·φ₀/Σw, exact
only at flat-flux equilibrium) → ``AngularEdgeExtrapolation``.
Post-fix the curvilinear isotropic MMS collapses O(h²) into the
absolute-magnitude band (sphere orders 2.00–2.01 to 5.74e-5 at
nx=320; cylinder 2.00 to 3.37e-5 at nx=160), SI ≡ Krylov
bit-identical, and the 2 isotropic xfail-strict tripwires in
``tests/sn/verification/mms/test_mms_curvilinear.py`` are REMOVED
(now plain tests carrying ``catches("ERR-058")``).  The 4 anisotropic
xfail markers were REMOVED by W3 (Issue #229, commit ``679a1e6``,
2026-06-13) — they covered a DIFFERENT, non-ERR-026 reason (the
fixed-quadrature angular half-angle-thread INTERPOLATION floor of the
per-ordinate-imposed aniso ansatz), and the gates were retuned to
assert what is TRUE (sphere coarse O(h²) window + magnitude band;
cylinder verified floor-scaling, no rate claim — the cylinder has no
pre-floor O(h²) window at any practical quadrature).  Phase A–D
narrative above preserved as history.  Closes #98 → #99 → #168 → #195.

**Surviving manifestation (distinct root, NOT the wrong-fixed-point
class).**  ERR-026 was the curvilinear *wrong-fixed-point* family
(closure-seed bias, terminally closed by ERR-058 / #195 / #196).  What
SURVIVES the curvilinear-aniso program (W1–W5, 2026-06-13) is a
DIFFERENT, inherent defect: the central-cell (r → 0) scalar flux is
first-order O(h) in the pointwise / L∞ norm because plain diamond is
inconsistent at the zero-area inner face — see **ERR-059** (Issue #233,
DOCUMENTED INHERENT LIMITATION, WONTFIX-for-DD).  It is NOT a residual
instance of ERR-026 (the fixed point is now CORRECT; the limitation is
the spatial *order* at the singular cell, invisible to the L2 / k_eff
gates that ERR-026 ever blocked).  The W1 ``τ``-clamp mis-citation
finding (sphere unclamp, commit ``b2d8a6d``) is the other W1–W5 leg of
this family's clean-up — see the ``τ-clamp mis-citation finding``
entry at the end of this catalog.

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
returns this constant. See ``docs/theory/references/peierls_nystrom.rst``
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

**Test reference (Wave 0 / Round 1):** `tests/numerics/test_projection_operators.py::TestApplyTransposeIsWWeightedAdjoint::test_adjoint_identity_matches_production` and `::test_apply_transpose_no_2l_plus_1_factor`. Tagged `@pytest.mark.l1` and `@pytest.mark.catches("ERR-039")`. *(2026-07-21 sync: this file and class no longer exist — the live catchers are the Round-2 reference below.)*

---

**2026-05-26 Phase 1 update (refactor/moment-space-and-layering):**

The Wave 0 fix made `apply_transpose` return the naked :math:`S_0(c)`
— but the Phase 1 audit revealed this was ALSO mislabeled. The bare
:math:`S_0` is **neither** the representation transpose **nor** the
Hilbert adjoint; both are :math:`S_0` post-multiplied by a different
diagonal:

- :math:`\Pi^\top = w_n \cdot S_0` (representation transpose) —
  `MomentProjection.apply_transpose` post-P1.4
- :math:`\Pi^* = g_C \cdot S_0` (Hilbert adjoint) — `MomentProjection.H`
  composed by the generic `_AdjointOperator` machinery using the
  codomain's SH-Gram metric
- :math:`R = (2\ell+1) \cdot S_0 = 4\pi \cdot g_C^{-1} \cdot S_0`
  (addition-theorem reconstruction) — `HarmonicMomentReconstruction`

All four operators (:math:`S_0, \Pi^\top, \Pi^*, R`) now have
separately-typed homes; conflating any pair is structurally
prevented. The :math:`(2\ell+1)` literal lives in exactly one place:
`SphericalHarmonicBasis.addition_theorem_factor`; the :math:`g_C =
\mathrm{diag}(4\pi/(2\ell+1))` metric on
`SphericalHarmonicSpace.inner_product_weights`.

The Wave 0 test `TestApplyTransposeIsWWeightedAdjoint` was renamed to
`TestApplyTransposeIsRepresentationTranspose` (still
`@pytest.mark.catches("ERR-039")`); a new class
`TestHilbertAdjointViaGenericMachinery` covers the :math:`g_C \cdot
S_0` identity. The new file
`tests/numerics/test_spherical_harmonic_space.py` (13 tests, all
`@pytest.mark.catches("ERR-039")` on the math-identity tests) carries
the Phase 1 endpoint.

**Lesson generalised:** When a bug is fixed by "removing the wrong
operator", verify that the replacement is the **right** operator under
the discipline's metric, not just "different from the wrong one". The
Wave 0 fix removed the (2ℓ+1) factor but did not install the W metric
that the W-weighted adjoint identity demands. The pattern was
"fix-the-symptom-not-the-cause" applied to a convention-drift root —
a deeper fix would have asked: "what IS the right adjoint?" rather
than "what was wrong about the old one?". See ERR-051 for the related
discipline failure on the validation-method side.

**Phase 1 test reference (Round 2; synced 2026-07-21):** `tests/numerics/test_spherical_harmonic_space.py` — the live catcher set (8 tests marked `@pytest.mark.catches("ERR-039")`). The adjoint endpoint is `test_H_equals_g_C_times_S0` (pins `M.H ≡ g_C · S₀`, the metric times the naked synthesis); `test_R_equals_2l_plus_1_times_S0` pins the complementary (2ℓ+1)-lives-in-R side. `tests/numerics/test_projection_operators.py` (with its `TestApplyTransposeIsRepresentationTranspose` / `TestHilbertAdjointViaGenericMachinery` classes) was retired tree-wide by the Frame campaign; its coverage lives in `test_spherical_harmonic_space.py`. Phase 1 commits: 0eb9cf3..c5be4b0 on `refactor/moment-space-and-layering`.


---

## ERR-040 — Tangential ordinate silently classified as inflow OR outflow at a face requiring strict partition

**Status:** **TYPE SHIPPED Wave 3 2026-05-11** (refactor/sn-operator-algebra
→ feature/bc-trace-law-abc). Concrete BC that fires the invariant lands
in Wave 7. This entry pins the error TYPE; the catching test ships
with the Wave 7 concrete that overrides
`assert_inflow_outflow_classification` to raise.

**Failure mode:** **#5 (index error)** — an ordinate with
`|Ω · n_f| <= eps` (a quadrature node that grazes the face within
floating-point tolerance) is inappropriately included in either the
inflow or the outflow half. The Wave 2 `_TANGENTIAL_EPS=1e-12`
already excludes such ordinates from BOTH masks in
`InflowTraceSpace` / `OutflowTraceSpace`, but a downstream consumer
that assumes the inflow + outflow partition is a strict cover of the
ordinate set will silently miscount the tangential ordinate.

**Module:** `orpheus.geometry.boundary._errors.IncomingOutgoingTraceClassificationError`.

**Mechanism:** Reflective BCs require a clean inflow→outflow
involution: every inflow ordinate at the face must have an outflow
partner under reflection. If a tangential ordinate (Lebedev
`(0, 0, ±1)` on the perpendicular face is the canonical case) is
treated as either inflow or outflow, the reflection partnership is
ill-defined — the partner is the tangential ordinate ITSELF. The
``IncomingOutgoingTraceClassificationError`` is the typed signal
that the face contract requires a strict cover and one is missing.

**How it hides:** A flat-flux test with uniform ψ across all
ordinates makes the tangential count irrelevant (its contribution
cancels). The ERR fires only under non-trivial angular dependence
AND a reflective BC contract. Caught by the Wave 7 invariant.

**Catching test (Wave 7):** the concrete `ReflectiveBoundary`
ships a test that constructs a Lebedev-13 quadrature against a
face whose normal aligns with the tangential ordinate (i.e.,
`mu_z` for a `zmin/zmax` face on a Lebedev quadrature without
the `(0, 0, ±1)` node guard), then verifies
`assert_inflow_outflow_classification` raises
`IncomingOutgoingTraceClassificationError`. To be tagged
`@pytest.mark.catches("ERR-040")` at Wave 7 ship.

**Lesson:** A strict-partition contract is part of the BC's TYPE
SIGNATURE. The trace-space layer (Wave 2) excludes tangential
ordinates from both masks; the law layer (Wave 3) must raise when
its contract requires a strict cover. The two layers' contracts
COMPOSE — neither alone is sufficient.


---

## ERR-041 — Vacuum BC constructed against an outgoing trace (Γ_+ instead of Γ_-)

**Status:** **TYPE SHIPPED Wave 3 2026-05-11.** Concrete realiser
ships Wave 5.

**Failure mode:** **#6 (convention drift)** — a vacuum BC defined as
"set γ_- ψ = 0 on the INCOMING trace" gets constructed against an
outgoing-trace metadata object (typically because a Wave 5 realiser
swapped the face's `(inflow, outflow)` annotation by mistake). The
resulting operator zeroes the outgoing flux — geometrically
meaningless and provably wrong.

**Module:** `orpheus.geometry.boundary._errors.VacuumAppliedToOutgoingTraceError`.

**Mechanism:** Vacuum is the rank-0 case of the affine form
`γ_- ψ = R G γ_+ ψ + q` with `R = 0`, `G = 0`, `q = 0`. The
operator's DOMAIN is the outgoing-trace space `Γ_+`; the RANGE is
the inflow-trace space `Γ_-`. Constructing the operator with
domain = `Γ_-` would make the apply path zero the inflow we just
computed.

**How it hides:** Bit-identical legacy tests (`apply` returns
`np.zeros_like(psi_out)` regardless of trace orientation) pass —
the bug only surfaces when the trace metadata is consumed by
Wave 5+ realiser logic that asserts `domain == outflow_trace`.

**Catching test — SHIPPED (coupled-block campaign task #52,
`51c22396`, 2026-07-12):** the invariant went PRODUCTION-live as the
realizer-seam guard in `SNBoundaryRealizer.realize`'s vacuum arm —
claimed `inflow_indices` cross-checked against the orientation the
FACE NAME alone implies via `build_omega_dot_n` (the canonical
`SNMesh.realize_boundary_law` path derives both from the trace ⟹
green by construction; the hand-built / annotation-swap path is
exactly where the two encodings are independently supplied).
Catchers tagged `@pytest.mark.catches("ERR-041")` in
`tests/sn/operators/test_sn_boundary_realizer.py` (swap /
single-index contamination / positive both faces / the documented
faceless escape / unknown-face loud), mutation-verified in-process
(no-op the guard ⟹ the raise disappears; revert ⟹ red).

**Lesson:** Trace orientation is part of the BC's type signature.
When the realiser layer (Wave 5) wires a BC to a trace, it MUST
assert the orientation matches the BC's contract — silent
misorientation produces a silent zero where flux belongs.


---

## ERR-042 — Reflection-index table inconsistent with quadrature weights (measure-non-preserving G)

**Status:** **TYPE SHIPPED Wave 3 2026-05-11.** Catching test
ships with `assert_geometry_map_measure_preserving` override in
Wave 7 ReflectiveBoundary.

**Failure mode:** **#5 (index error)** + **#6 (convention drift)**
— a reflection table built against `mu_x` but applied to a
quadrature whose ordering treats `mu_x = -mu_y` (or any other
permutation drift) fails to preserve the direction-cosine measure
`w(Ω) · |Ω · n|`.

**Module:** `orpheus.geometry.boundary._errors.BoundaryGeometryMapNotMeasurePreservingError`.

**Mechanism:** A reflection that does not preserve the measure
means the reflected current at the face differs from the incident
current — physically impossible for an isometry (specular reflection
is an isometry). The bug shows up as a small but persistent net
current at a face where there should be none (zero leakage from a
purely reflective boundary).

**How it hides:** A homogeneous problem with uniform Σ_t hides
the measure mismatch because the imbalance is small (a few ULP per
ordinate times the weight). Heterogeneous problems amplify the
imbalance by the contrast in Σ_t.

**Catching test (Wave 7):** the ReflectiveBoundary's
`assert_geometry_map_measure_preserving` got its REAL body at the
coupled-block campaign task #52 (`51c22396`, 2026-07-12): the trace
measure `m = w·|μ_axis|` must satisfy `m[π] ≈ m` (rtol 1e-12),
INDEPENDENT of involution — the pre-#52 delegation ("weight
equality implied by construction") left the
involutive-unequal-weight-pairing hole open, and the GL-8
neighbor-pair mutant passes involution while redding measure
(independence measured, both directions). Fired at realization by
`BoundaryTraceLaw.assert_realizable` (every `SNMesh` construction
certifies). Catchers tagged `@pytest.mark.catches("ERR-042")` in
`tests/geometry/test_bc_universal_invariants.py`
(positive+negative-paired, mutation-verified in-process).

**Lesson:** Geometric maps `G` carry a measure-preservation
contract. The contract must be CHECKED at construction, not only
verified by downstream global-current balance — by then the
miscount has propagated.


---

## ERR-043 — Boundary response kernel produces negative output (e.g. negative albedo via sign-flipped construction)

**Status:** **TYPE SHIPPED Wave 3 2026-05-11.** Catching test ships
with Wave 7 White / Albedo concretes.

**Failure mode:** **#1 (sign flip)** — the response kernel `R` is
constructed with a sign error (e.g. `-1.0 * G_diff` because of a
wrong unit-conversion factor), producing negative reflected current
from positive incident current.

**Module:** `orpheus.geometry.boundary._errors.BoundaryResponseNotPositiveError`.

**Mechanism:** Physically meaningful BCs (white, albedo) have
non-negative response kernels — the reflected flux is a fraction
(0 ≤ α ≤ 1) of the incident, never opposite-signed. A negative R
produces negative ψ at the inflow trace, which propagates through
the sweep as locally negative flux. The bug masquerades as
"negative-flux artifact"; root cause is the sign of R.

**How it hides:** Eigenvalue solves can mask sign errors by
normalising the dominant eigenvector to be non-negative — the
sweep's local negative-flux artifact looks like a "discretisation
artifact" and is silently fixed.

**Catching test (Wave 7):** White / Albedo BC tests that probe
the response on a known positive incident flux and assert the
output is element-wise non-negative; the `assert_response_positive_if_declared`
invariant fires the error. To be tagged
`@pytest.mark.catches("ERR-043")` at Wave 7 ship.

**Lesson:** Sign invariants on response kernels MUST be asserted at
construction. Downstream consumers (sweeps, Krylov) cannot
distinguish a wrong-sign R from a discretisation artifact.


---

## ERR-044 — Reflection permutation is not an involution (perm ∘ perm ≠ id)

**Status:** **TYPE SHIPPED Wave 3 2026-05-11.** Catching test
ships with Wave 7 ReflectiveBoundary invariant override.

**Failure mode:** **#5 (index error)** — the reflection-index
construction produces a permutation that, applied twice, does NOT
return every ordinate to itself. Typically a row-order vs
column-order bug or a partial table (some ordinates missing from
the reflection map).

**Module:** `orpheus.geometry.boundary._errors.ReflectionNotInvolutiveError`.

**Mechanism:** An axis reflection is its own inverse by
definition (apply twice → identity). The
`PermutationOperator.is_involution` flag (Wave 0) detects this at
construction. If the involution check fails, the BC's
`apply_transpose` returns a different permutation than `apply`,
breaking the adjoint identity that sensitivity-analysis pipelines
depend on.

**How it hides:** Forward-only solves never exercise
`apply_transpose`. A non-involutive reflection table passes every
forward eigenvalue test and only fails when an adjoint pipeline
(generalised perturbation theory, sensitivity coefficients) is
introduced.

**Catching test (Wave 7):** ReflectiveBoundary test that
constructs a deliberately-broken reflection table (drop a row,
swap two indices), passes it to the BC, and asserts
`ReflectionNotInvolutiveError` fires from
`assert_geometry_map_measure_preserving` or a dedicated involution
assertion. To be tagged `@pytest.mark.catches("ERR-044")` at
Wave 7 ship.

**Lesson:** Involution is a constructional invariant on the
reflection table. The Wave 0 `PermutationOperator.is_involution`
flag exists for this purpose; the BC layer must consult it and
raise on failure rather than silently allowing a non-involutive
permutation.


---

## ERR-045 — Reflection maps an inflow ordinate to itself rather than to its outflow partner

**Status:** **TYPE SHIPPED Wave 3 2026-05-11.** Catching test
ships with Wave 7 ReflectiveBoundary invariant override.

**Failure mode:** **#5 (index error)** — for a non-axis-aligned
reflection (or an axis-aligned reflection on a quadrature node
near the axis), an inflow ordinate at the face is mapped back to
itself, not to its outflow partner.

**Module:** `orpheus.geometry.boundary._errors.ReflectionDidNotMapInflowToOutflowError`.

**Mechanism:** Specular reflection's geometric contract is
"every inflow Ω_n at the face has a partner Ω_n' with
Ω_n · n = -Ω_n' · n". A reflection that maps inflow to inflow
sets up a self-loop in the sweep dependency graph that the
`creates_sweep_cycle` flag was meant to flag — but the cycle is
TRIVIAL (length 1) and degenerates the sweep convergence theory.

**How it hides:** Self-loops at a tangential ordinate (or a
nearly-tangential one) look like a benign quadrature artifact and
don't surface until per-ordinate flat-flux residual tests stress
the redistribution path (cf. ERR-006 / Signature 1).

**Catching test — SHIPPED (coupled-block campaign task #52,
`51c22396`, 2026-07-12):** the new
`assert_reflection_maps_inflow_to_outflow` invariant — every
non-tangential ordinate's partner must be non-tangential with
OPPOSITE sign on the law's axis (tangential exemption rides the one
shared `TANGENTIAL_EPS`, promoted public); raises
`ReflectionDidNotMapInflowToOutflowError`. The identity-table
mutant passes involution AND measure and reds ONLY here — the
three-independent-invariants lesson realized as three MEASURED
independent checks, fired at realization via
`BoundaryTraceLaw.assert_realizable`. Catchers tagged
`@pytest.mark.catches("ERR-045")` in
`tests/geometry/test_bc_universal_invariants.py`
(positive+negative-paired incl. realize-time legs poisoning
`quad.reflection_partners`, mutation-verified in-process).

**Lesson:** Reflection contracts compose: the inflow partition,
the involution property, and the inflow → outflow image are three
independent invariants. All three must hold; checking only one or
two leaves a hole.


---

## ERR-046 — Albedo / white kernel with α > 1 (sub-Markov violation)

**Status:** **TYPE SHIPPED Wave 3 2026-05-11.** Catching test
ships with Wave 7 Albedo / White concretes.

**Failure mode:** **#4 (factor error)** — an albedo value
greater than 1 is admitted by the constructor. A row-sum of the
response kernel `R` strictly exceeds 1, violating the sub-Markov
property that physically describes a non-amplifying boundary.

**Module:** `orpheus.geometry.boundary._errors.SubmarkovViolationError`.

**Mechanism:** `∫ R dy ≤ 1` is the canonical sub-Markov condition
on a stochastic kernel (cf. Bell & Glasstone 1970 §1.5). An
albedo greater than 1 implies a source on the boundary surface,
which the BC framework does NOT model (sources live in `q`, not
in `R`). Admitting α > 1 can drive the eigenvalue iteration to a
spurious supercritical state where the source comes from the
boundary, not the volumetric fission.

**How it hides:** k-eigenvalue iteration rescales the flux to
the dominant eigenvector regardless of where the "source" lives;
the wrong-α run produces a higher k and a flux profile peaked
near the boundary — both look physically plausible to anyone not
checking the input data.

**Catching test (Wave 7):** Albedo constructor test that passes
`α=1.2` and asserts `SubmarkovViolationError` fires. White
constructor test that uses a kernel construction known to row-sum
to 1.5 and asserts the same error. To be tagged
`@pytest.mark.catches("ERR-046")` at Wave 7 ship.

**Lesson:** Range constraints on physical parameters MUST be
enforced at construction. "α is a float" is not a contract; "α ∈
[0, 1]" is. The typed error makes the contract grep-able and
testable.


---

## ERR-047 — Boundary source q has nonzero entries on the outflow trace (q ∉ Γ_-)

**Status:** **TYPE SHIPPED Wave 3 2026-05-11.** Catching test
ships with Wave 7 prescribed-inflow BC.

**Failure mode:** **#6 (convention drift)** — the source `q` in
the affine form `γ_- ψ = R G γ_+ ψ + q` is constructed with
nonzero entries on outflow-trace ordinates. The downstream
realiser then writes a positive flux into outflow slots that the
sweep will OVERWRITE — silently dropping the source contribution.

**Module:** `orpheus.geometry.boundary._errors.BoundarySourceNotOnIncomingTraceError`.

**Mechanism:** The affine form's `q` lives on `Γ_-` by definition.
If the user passes an array shaped like the full ordinate set but
populated uniformly (e.g. via `ConstantInflowSource(value=2.0)`
without masking), the outflow ordinates get spurious `+2.0`
entries that the sweep then discards. The total inflow current is
SHORT by the masked-out fraction — a missing-source bug that
masquerades as "boundary source is weaker than expected".

**How it hides:** A test that imposes a uniform inflow on EVERY
ordinate (including outflow) sums to the right total only because
the sweep silently zeroes the outflow side. A finer test that
sums per-ordinate inflow against the expected `q(Ω)` profile
catches it.

**Catching test (Wave 7):** prescribed-inflow BC test that
constructs `ConstantInflowSource(2.0)` on a quadrature with
mixed inflow / outflow ordinates: SHIPPED at the coupled-block
campaign task #52 (`51c22396`, 2026-07-12) — the ABC no-op default
became the REAL universal body (probe `self.source` on `(N,)`;
q ≡ 0 certifies trivially; nonzero q WITHOUT an inflow set raises
— "the realiser asked to apply the source without an outflow
mask", this entry's catcher spec verbatim), and the realizer's
`PrescribedInflow` arm now delivers q MASKED to Γ₋
(`IncomingSourceOperator` gained the inflow-indices plumbing;
delivered q = the source value on inflow slots, exactly 0 on
outflow AND tangential slots, pinned end-to-end — the three
pre-#52 tests blessing the unmasked full-array delivery, the
ERR-047 hazard itself, rewired to the masked contract). Catchers
tagged `@pytest.mark.catches("ERR-047")` in
`tests/geometry/test_bc_universal_invariants.py`
(positive+negative-paired, mutation-verified in-process).

**Lesson:** The affine form's three terms (`G`, `R`, `q`) each
live on a SPECIFIC trace space. The trace-space type tag (Wave 2)
+ the typed error (Wave 3) together encode this; future Waves
must respect the trace-space contracts rather than treating
`q` as "just a float vector".

---

## ERR-048 — Curvilinear SI sweep: pole-face WDD IC + Carlson seed source convention drift between SI and apply-matvec paths

**Date:** 2026-05-13 (Issue #196 Phase G Step 2 Path C).

**Status:** **CLOSED.** Two surgical patches landed in
`orpheus/sn/sweep.py::_sweep_1d_spherical` and `_sweep_1d_cylindrical`;
catching test `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`
(L0 streaming-equilibrium gauntlet, sphere + cylinder) ships
green; 5 curvilinear regression snapshots regenerated under the
corrected SI with three-pillar attestation; Phase E flux-shape
sentinel xfail-strict marker removed.

**Failure mode:** **#6 (convention drift)** — TWO simultaneous
convention disagreements between the SI sweep
(`_sweep_1d_spherical/cylindrical`) and the apply-matvec
(`transport_operator_matvec_spherical/cylindrical`).  Both
production paths converged to *their own* fixed point to machine
precision; the two fixed points differed by ~22% L0 on the
homogeneous reflective sphere streaming-equilibrium test.

**Module:** `orpheus.sn.sweep`.

**Mechanism:**

1. **Pole-face WDD IC (DOMINANT, ~95% of the L0 error).**
   For outward sweeps (μ ≥ 0) starting at the pole/axis cell
   i=0, the SI sweep initialised `ψ_face_in = 0` (`sweep.py:559`
   pre-fix).  The apply-matvec uses `ψ_face_in = ψ_cell[i=0]`
   per Lewis-Miller §4.5 Carlson starting-direction convention
   (`operator.py:781-786`).  At the pole face A[0] = 4π·0² = 0,
   so the streaming-IN term `|μ|·A[0]·ψ_face_in` is identically
   zero either choice — but the WDD diamond recurrence
   `ψ_face_out = 2·ψ_cell − ψ_face_in` propagates the choice
   downstream.  Setting `ψ_face_in = ψ_cell` preserves
   `ψ_face_out = ψ_cell` on a flat ψ field, hence Pomraning
   isotropy at the pole.  The SI choice `ψ_face_in = 0` produces
   `ψ_face_out = 2·ψ_cell` which destroys flat-field invariance
   from cell 0 onward.

2. **Carlson seed `Q_bar` source (sub-dominant, ~3% rel residual
   after fix 1 alone).**  The SI sweep built the Carlson coupled-
   pole seed source as `Q̄ = 0.5 · Q_within_group`
   (`sweep.py:514` pre-fix, where `Q_within_group` is scattering
   + fission/k + external).  The apply-matvec builds
   `Q̄ = 0.5 · Σ_t · φ_0(ψ_input)` (`psi_half_angle_seed.py:569`).
   At the within-group fixed point `Σ_t · φ_0 = Q_within_group`
   exactly, so the two conventions AGREE there; off the fixed
   point during Picard iteration they DIFFER.  Driving the SI's
   Carlson seed from the previous-iter scalar flux makes the SI
   Picard recursion consume the same seed as the apply-matvec.

**Both fixes carry previous-Picard-iter state through the
`psi_bc` dict** (existing infrastructure):

- `psi_bc["psi_pole"]` (sphere) / `psi_bc["psi_pole_cyl"]`
  (cylinder): `(N, ng)` cell-centre ψ at i=0 from the previous
  outer iteration; carries the pole-face IC.
- `psi_bc["phi_0_prev"]` (sphere) / `psi_bc["phi_0_prev_cyl"]`
  (cylinder): `(nx, ng)` scalar flux from the previous outer
  iteration; carries the Carlson seed source.

**Cold-start safety net:** on the first outer iteration neither
key exists, so both fixes degrade gracefully to the legacy
Phase F values (zero pole IC, `Q̄ = 0.5 · Q_within_group`).
Subsequent iterations use the apply-matvec conventions and
Picard converges to the apply-matvec fixed point.

**How it hid:**

- **Phase F's diagnostic compared SI-vs-Krylov k_eff** (agreed to
  0.286% at n=40, decaying as O(h)); the absolute magnitude of
  the convergence gap looked like "discretisation error
  approaching zero" rather than a convention bug.
- **Phase F regression snapshots were generated by the buggy SI**
  (`_generate_snapshots.py` uses default `inner_solver`); bit-
  identity preserved the wrong fixed point as the "regression
  truth" for 6 curvilinear snapshots.
- **L0 streaming-equilibrium was not tested**; the canonical
  Pomraning isotropy test `ψ_n(r=0) = Q/(Σ_t·Σw·(1−c))` was
  missing from the SN suite for the SI/sweep path.  The L1
  Gate 1.1 sphere MMS test had been adapted in Phase D for an
  isotropic-flat trial that masks the redistribution path
  (`vv-principles` §"Mode 7 — MMS simplification bias").
- **The legacy SI SOURCE iterates** (the Picard loop is on the
  scattering source, not the angular flux), so the Picard
  residual norm reaches machine precision on the wrong fixed
  point with no observable divergence.

**Which test catches it:**
`tests/sn/spatial/test_streaming_equilibrium_curvilinear.py` —
26 parametrised cases (sphere/cylinder × {20,40,80} cells ×
{4,8 or 8,16} ordinates × `{source_iteration, krylov}`) +
Pomraning pole isotropy gate (cv < 0.01).  Tagged
`@pytest.mark.l0 @pytest.mark.catches("ERR-048")
@pytest.mark.verifies("hebert-3-432")`.

**Lesson:** In curvilinear SN with pole/axis geometry, the
pole-face WDD initial condition for outward sweeps at i=0 MUST
mirror the cell-centre value (Lewis-Miller §4.5 Carlson
starting-direction); the Carlson coupled-pole seed `Q̄` MUST be
built as `0.5 · Σ_t · φ_0(ψ_prev)`, not `0.5 · Q_within_group`.
The apply-matvec used the correct conventions since Phase D
(Issue #168); the SI sweep used different conventions for both,
and the discrepancy manifested as mesh-independent 22% L0 pole
error.

**This is also a Pattern 2 (single source of truth) anti-pattern
instance:** the SI sweep and apply-matvec are TWO implementations
of the SAME continuous operator.  Twin paths drift; Phase D's
Carlson backport into the apply-matvec produced the
manifestation #6 pair that Phase F partially closed and Phase G
Step 2 surgically completes.  The architectural fix — one
shared per-cell operator (`SNCellOperator`) composed in both
call sites — lands in Phase G Steps 3+.

→ `numerical-bug-signatures` Signature 1 (curvilinear sweep
divergence under refinement) catalogues the manifestation
fingerprint family; ERR-026 (the original Phase D fix on the
apply-matvec) is this defect's twin.

### Manifestation #3 — cylinder per-level Carlson seed Q_bar normalisation hardcoded for sphere Σw=2

**Date:** 2026-05-13 (Issue #196 Phase G Step 2 cylinder fix).

**Failure mode:** **#6 (convention drift)** — sphere quadrature
normalisation `Σw_sphere = 2` was hardcoded as the literal `0.5`
into the Carlson seed kernel ``Q_bar = 0.5 · Σ_t · φ_0``.  For
sphere this literal happens to equal ``1/Σw_full`` and produces
the correct seed; for any 3D quadrature (ProductQuadrature,
LevelSymmetric, Lebedev — all with ``Σw_full = 4π``) the literal
is wrong by factor ``2π``, producing a divergent residual.

**Module:** `orpheus.sn.spatial.psi_half_angle_seed` (apply-matvec
path), `orpheus.sn.sweep` (SI sphere + cylinder paths).

**Mechanism:** The Hébert §3.9.4 (3.432)-(3.435) coupled-pole
inward μ = −1 sweep computes the half-angle seed from a cell-
averaged source ``Q̄ = (2ℓ+1)/2 · Σ_t · φ_0 · P_0(−1)`` for the
``ℓ = 0`` (isotropic) Legendre moment.  The factor ``(2ℓ+1)/2``
is the sphere-1D angular measure ``∫P_0² dμ / 2 = 1/2`` valid
on ``dμ`` over ``[−1, 1]``.  For 3D quadratures the
normalisation that preserves Pomraning isotropy is
``1/Σw_level`` (apply-matvec, where ``φ_0`` is the per-level
moment) or ``1/Σw_full`` (SI sweep, where ``phi_0_prev`` is the
full scalar flux); for sphere these coincide at ``1/2`` and the
sphere-only literal is bit-identical to either.  For cylinder
ProductQuadrature the per-level azimuthal weight sum is
``2π`` ≠ ``2``, and the literal ``0.5`` produces a seed
overshooting by ``2π × 0.5 = π``.

**Empirical pre-fix evidence:** ``L · ψ_flat − Σ_t · ψ_flat`` on
homogeneous reflective cylinder mixture B 1G:

```
n_mu=2, n_phi=2:  n_cells=20: ||res|| = 6.03   → n_cells=160: 17.2 (GROWS)
n_mu=4, n_phi=4:  n_cells=20: ||res|| = 3.40   → n_cells=160: 9.71 (GROWS)
```

End-to-end fixed-source error: 580 % at ``[source_iteration-4-20]``
gauntlet case.  Signature 1 (curvilinear sweep divergence under
refinement) fingerprint.

**How it hid:**

- The Phase D Carlson seed was developed on sphere first; the
  literal ``0.5`` was correct for sphere's GL-N quadrature and
  the implementation was carried forward to cylinder unchanged.
- The cylinder per-level loop introduced
  ``CarlsonSweepContext.weights = weights_level`` (per-level
  weight sub-array) — but ``CarlsonInwardSweep.__call__`` kept
  the sphere literal instead of consuming
  ``1/context.weights.sum()``.  The defect lived on a single
  shared file:line (``psi_half_angle_seed.py:569``).
- Phase G Step 2 Path C **claimed** "12/12 cylinder PASS"
  in its closeout memo without running the cylinder test; the
  empirical run (this manifestation's catching test) shows
  12/12 FAILURE with up to 580 % rel error at the smallest case.

**Empirical post-fix evidence:**

```
Cylinder L0 streaming-equilibrium gauntlet (12 cases):
  source_iteration, n_mu=4/8, n_cells=20/40/80:  all PASS at rtol=1e-9
  krylov,           n_mu=4/8, n_cells=20/40/80:  all PASS at rtol=1e-9

Convergence under refinement (post-fix, production solver):
  n_cells=20→160, n_mu=2/4, n_phi=2/4:  err ~ 1.5e-11 (mesh-INDEPENDENT)
  Krylov ≡ SI ≡ analytical to ~2e-12 on every grid

Sphere streaming-equilibrium gauntlet (12 cases):  all PASS  (no regression)
Apply-matvec invariant ``L·ψ_flat = Σ_t·ψ_flat`` (12 cases):  PASS at FP-noise
Phase E flux-shape sentinel (2 heterogeneous MR cases):  PASS
SN regression suite (11 snapshots, including regen'd cylinder):  all PASS
```

**Which test catches it:**

- `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py::test_homogeneous_streaming_equilibrium_cylinder`
  — the canonical L0 streaming-equilibrium gauntlet (12 cases).
- `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py`
  — NEW.  Promoted from the numerics-investigator's diagnostic
  ``derivations/diagnostics/diag_phase_g_step2_cyl_residual_pytest.py``.
  Covers (a) the apply-matvec flat-flux invariant (12 cases) +
  (b) the 4-leg 3-way standoff (12 cases): Krylov vs analytical,
  SI vs analytical, Krylov ≡ SI at machine precision.

**Defect sites (now fixed):**

- `orpheus/sn/spatial/psi_half_angle_seed.py:569` (apply-matvec
  Carlson seed kernel): ``0.5 · Σ_t · φ_0`` →
  ``Σ_t · φ_0 / weights.sum()``.
- `orpheus/sn/sweep.py:543/547` (sphere SI sweep): ``0.5 · …`` →
  ``… / weights.sum()`` (bit-identical for sphere; hygiene fix
  per ``coding-elegance`` Pattern 7).
- `orpheus/sn/sweep.py:754/756` (cylinder SI sweep): same form.

**Three-pillar attestation for regenerated cylinder snapshot
(``cyl_1g_homogeneous_product_dd_n20``):**

1. **L0 streaming-equilibrium**: regenerated keff = 1.5 (exact
   k_∞ = νΣ_f/Σ_a for mixture A 1G), scalar_flux = 1.5 flat
   (std = 1.6e-11, cv = 1.0e-11) — Pomraning isotropy + reflective
   symmetry on homogeneous infinite medium.
2. **Pomraning pole isotropy**: cv(ψ@i=0) ≪ 0.01 (gate by ~9
   orders of magnitude).
3. **Variant α cross-check (structurally independent reference)**:
   ``test_phase_e_trajectory_resolvent_flux_shape_crosscheck`` 2/2
   PASS at Phase E rtols (``cyl_2g_3reg_LS4_dd_n40`` heterogeneous
   MR cylinder geometry passes — the new SI fixed point matches
   the bouncing-characteristic Variant α reference shape).

**Lesson:** Numerical literals that LOOK like ``0.5`` may be
encoding a quadrature-dependent normalisation
``1/Σw_quadrature``.  Replace all such literals with
``1.0/weights.sum()`` from the quadrature-in-scope per
``coding-elegance`` Pattern 7 (normalise at the definition site,
not at every consumer).  The literal worked for sphere by
coincidence (Σw=2) and broke silently on cylinder
(Σw_level=2π).  This is a Pattern 7 failure that the original
Phase D implementation should have avoided by construction.

**Operating principle reaffirmed:** verification claims **MUST**
be backed by VERBATIM paste-back from a real test run.  The
prior Phase G Step 2 Path C closeout's "12/12 cylinder PASS"
claim was unverified; the cylinder test was never actually run
during that step.  Hallucinated test results are a session
failure per Cardinal Rule 1 (correctness is CRITICAL).

→ `numerical-bug-signatures` Signature 1 (curvilinear sweep
divergence under refinement) catalogues this fingerprint;
ERR-026 manifestations #6 + #7 are this defect's twins on
sphere (closed by Phase F+Phase G Step 2 Path C).


---

## ERR-049 — Convention drift between operator-algebra and transport_sweep — per-ordinate vs iso magnitude on the typed AngularFlux carrier

**Date:** 2026-05-19 (R-1 Step 4 session 1 Step E — SI carve onto
`SourceIteration` + typed `AngularFlux` diagnostic).

**Status:** **CLOSED.**  Architectural fix landed in Phase 1.1 / A1
commit ``de8822d`` (2026-05-22) on
``refactor/sn-operator-algebra``: producer-side ``/sum_w``
normalisation at every site that emits a per-ordinate value
(``ScatteringOperator.apply`` AngularFlux variant,
``FissionOperator.apply`` AngularFlux variant,
``PerOrdinateSource.from_isotropic`` factory).  The legacy ``×W``
bridge in ``StreamingCollisionOperator.solve`` and the ``/W`` rescales in
``_solve_krylov`` / ``_solve_source_iteration`` /
``_make_sweep_preconditioner`` all dissolved.

**Failure mode:** **#6 (convention drift)** — the operator-algebra
layer (``StreamingCollisionOperator``, ``ScatteringOperator``,
``FissionOperator``) carried values in **per-ordinate units**
(``ψ_n``, ``q_n = q_iso/Σw``); ``transport_sweep`` consumed an
``aniso_source`` argument in **iso magnitude** and applied an internal
``weight_norm = 1/W`` (``sweep.py:432`` pre-fix).  The carve in R-1
Step E threaded per-ordinate ``rhs.values`` directly through
``StreamingCollisionOperator.solve`` to ``transport_sweep``: the per-ord
``q_n`` got divided by ``W`` a second time, planting the converged
fixed point at ``k_inf/W`` instead of ``k_inf``.

**Module:** `orpheus.sn.operator.StreamingCollisionOperator.solve` (the
consumer-side bridge site); root cause is the convention asymmetry
between `orpheus.sn.scattering.ScatteringOperator.apply` (typed
producer; pre-A1 emitted iso magnitude on the AngularFlux variant)
and `orpheus.sn.sweep.transport_sweep` (consumer; pre-A1 expected
iso magnitude and applied internal ``/W``).

**Mechanism:** the operator-algebra factorisation
``(L + C − S)·ψ = q_ext`` is unit-stable at the typed level — every
operand carries the per-ordinate magnitude ``ψ_n`` and the per-
ordinate source ``q_n = q_iso/Σw``.  The pre-A1 ``transport_sweep``
took the iso-magnitude scalar source ``q_iso`` and applied
``weight_norm = 1/Σw`` internally as part of the within-cell balance
``q_iso/Σw + Σ_t · ψ_n``.  The two conventions are individually
correct; the bug is at the **boundary** where the operator-algebra
hands its per-ordinate source to the sweep, expecting "this is
already per-ord" and the sweep applying ``/Σw`` anyway.

For a homogeneous reflective uniform-source problem
``Σ_t · ψ_n = q_n``, the *correct* fixed point is
``ψ_n = q_n/Σ_t = (q_iso/Σw)/Σ_t``.  Pre-fix the SI carve produced
``ψ_n = q_n/(Σw · Σ_t) = q_iso/(Σw² · Σ_t)`` — a factor of ``Σw``
too small.  For an eigenvalue problem the same factor propagates to
``k_eff = (νΣ_f·φ) / (Σ_a·φ)`` indirectly via the within-group SI
recursion's source magnitude; the empirical drift was
``k_recovered/k_ref ≈ 1/(Σw·c_s)`` where ``c_s = Σ_s/Σ_t`` is the
group's scattering ratio.  For Gauss-Legendre N=8 (``Σw = 2``) and a
typical scattering-dominated thermal group ``c_s ≈ 0.9``, the
observed ratio was ``1/(2·0.9) = 0.556`` — consistent with the
``1.4844/1.875 = 0.7917`` slab-2eg signature once the eigenvalue
power-iteration's group-source mixing is accounted for.

**How it hid:**

- **Krylov inner solver was unaffected.**  The Krylov path
  (``_solve_krylov``) called the sweep ONLY as a preconditioner
  (``preconditioner=lambda q: q`` per session 1 Step D's explicit
  identity), so the convention bridge was a no-op.  The eigenvalue
  outer iteration converged for the Krylov path on every
  ``homogeneous_kinf`` case at rtol ≤ 1e-8.  Slab-2eg-SI was the
  failing leg.

- **Slab-1eg-SI passed by accident.**  One-group problems are
  flux-shape independent (``k_inf = νΣ_f/Σ_a`` doesn't depend on
  ``φ`` shape per the ``vv-principles`` 1-group-degeneracy rule).  The
  convention drift left the relative shape of ``φ`` unchanged on
  homogeneous-1g; the keff ratio between corrected and drifted
  paths cancels.  Multi-group eigenvalue is where the bug surfaces.

- **L0 streaming-equilibrium tests were missing for the SI path**
  in the pre-Phase-1.1 suite.  The L1 ``test_kinf_homogeneous.py``
  was the only multi-group gate; its first failing case
  (``slab-2eg-source_iteration``, ``keff=1.4844`` vs
  ``ref=1.875``) was the diagnostic trigger.

- **Pre-A1 the consumer-side ``×W`` bridge in
  ``StreamingCollisionOperator.solve`` had been added in R-1 Step C to
  paper over the asymmetry**; Step E's principled carve removed
  the bridge and exposed the underlying convention drift.

**Empirical pre-fix evidence (R-1 Step E pre-Phase-1.1):**

```
Slab GL-N=8 (Σw = 2.0):
  source_iteration:  keff = 1.4844  ref = 1.8750  rel_err ≈ 21%
  krylov:            keff = 1.8750  ref = 1.8750  rel_err ≈ 0% (unaffected)

Sphere GL-N=8:
  source_iteration:  keff = 1.4844  ref = 1.8750  rel_err ≈ 21%

Cylinder ProductQuadrature (Σw_full = 8.0):
  source_iteration:  keff scaled by similar factor

Streaming-equilibrium fixed-source on slab:
  pre-fix:   ψ_n = q_iso / (W² · Σ_t)   (factor W too small)
  post-fix:  ψ_n = q_iso / (W · Σ_t)    (correct)
```

**Empirical post-fix evidence (Phase 1.1 A1):**

```
L1 kinf_homogeneous suite (8 coords × 2 inner_solvers × 3 ng_keys + 6 cyl):
  35 passed + 5 xfailed (xfails pre-existing per #200)

L1 bridge-regression suite (Phase 1.3, tests/sn/test_invertible_operator.py
TestInvertibleSolveBridgeRegression):
  10/10 PASS — including all 3 coords × {2eg, 4eg} SI carve k_inf at rtol < 1e-9.

L0 streaming-equilibrium curvilinear:
  77 pass + 1 skip on the spatial/ slow gate
```

**Which test catches it:**

`tests/sn/test_invertible_operator.py::TestInvertibleSolveBridgeRegression`
(R-1 Step 4 Phase 1.3 promotion from
``derivations/diagnostics/diag_r1_step_e_invertible_solve_w_bridge.py``).
Three test methods, ``@pytest.mark.l1 @pytest.mark.catches("ERR-049")
@pytest.mark.verifies("transport-cartesian",
"sn-curvilinear-homogeneous-kinf-recovery")``:

1. ``test_slab_uniform_roundtrip`` — slab
   ``(L+C).solve((L+C).apply(ψ=1)) == ψ=1`` to machine zero
   (streaming-free check; pre-fix gave ``1/W``).
2. ``test_fixed_source_homogeneous_reflective_recovers_q_over_sigma``
   — parametrised over slab/sphere/cylinder; converges to
   ``ψ_n = q_n/Σ_t`` per ordinate (pre-fix gave ``q_n/(W·Σ_t)``).
3. ``test_si_carve_recovers_analytical_kinf`` — parametrised over
   {slab, sphere, cylinder} × {2eg, 4eg}; end-to-end SI inner
   solver on homogeneous reflective converges to analytical
   ``k_inf`` at rtol < 1e-9 (pre-fix failed every ng≥2 case at
   ~21% drift).

The N5 unit pin
``tests/sn/test_invertible_operator.py::TestSolve::test_solve_consumes_per_ordinate_rhs``
(``@pytest.mark.l0``) structurally pins the new contract:
``StreamingCollisionOperator.solve`` forwards ``rhs.values`` to
``transport_sweep`` bit-equal — no ``×W`` or ``/W`` rescaling on the
hot path.

**Lesson** (cf. `.claude/lessons.md` L17 + L18): a carve that crosses
subsystem boundaries (operator algebra ↔ sweep, scalar ↔ per-ord,
packed ↔ typed) MUST start with a written convention crosswalk
table.  When the producer and consumer disagree, fix the
**producer** (Pattern 7 per ``coding-elegance``): place the
``/Σw`` normalisation at the producer boundary, not at every
consumer.  Consumer-side bridges multiply with every new consumer;
producer-side normalisation costs once and dissolves the bridge
class entirely.

**Why a capability-style fix wouldn't have worked:** during R-1
session 1 a brief flirtation considered "the convention is part of
the value's type — wrap ``q`` in a ``PerOrdinateSource`` vs
``IsotropicSource`` and let dispatch decide".  The implementation
sketch ran into the dispatch-on-magnitude trap: the underlying ψ
values are still floats, so the type tag only documents the
expected magnitude, doesn't enforce it.  The principled fix is
to normalise at the producer — once at the boundary, never again
downstream.

→ `numerical-bug-signatures` Signature TBD (multi-group SI keff
drift by factor ``1/(Σw·c_s)`` while Krylov passes) catalogues
the fingerprint family; ERR-048 (curvilinear SI sweep convention
drift) is the structural sibling on a different convention axis
(Carlson seed source vs WDD pole-face IC, not per-ord magnitude).


---

## ERR-050 — Silent preconditioner fallback breaks stateful-inverse contract — `KrylovAcceleration(preconditioner=None)` defaulted to `L.solve` which read history from GMRES residual vectors

**Date:** 2026-05-19 (R-1 Step 4 session 1 Step D — Probe B sphere
identity-precond diagnostic).

**Status:** **CLOSED via structural supersession.**  The original
session-1 plan A2 proposed advertising a ``CAP_STATELESS_INVERSE``
capability and gating ``KrylovAcceleration`` on it (issue #203).
On 2026-05-22 the user-led architectural read identified the
deeper cause: 1D ``transport_sweep`` was duplicating the M-M
Carlson coupled-pole seed inline (calling
``carlson_inward_sweep_from_source`` directly + computing
``Q_bar`` inline) while the matvec already routed through
``MorelMontryAngularSweep.psi_half_seed``.  Phase 1.2
(commit ``c93355c``) unified both consumers through M-M's
``psi_half_seed`` strategy: ``StreamingCollisionOperator.solve(rhs, *,
initial_guess=None)`` is now a pure function — the Carlson seed
travels through the explicit ``initial_guess`` kwarg, not through
the lag-1 frame ``rhs(1)``.  ``KrylovAcceleration``'s
``preconditioner=None`` default invokes ``L.solve(q)`` with
``initial_guess=None``, which the M-M closure interprets as
**explicit cold start** (deterministic zero seed via
``psi_half_seed`` cold-path).  The silent-fallback path no longer
exists; the capability-flag patch is unnecessary.  Issue #203
closed by supersession.

**Failure mode:** **#4 (wrong recursion / state)** — a "stateless"
fallback path silently coupled to caller state that GMRES residual
inputs did not carry.  The bug class is "the default invokes a
primitive whose precondition (stateful caller) is not advertised
in the type system".

**Module:** `orpheus.sn.operator.StreamingCollisionOperator.solve` (pre-fix
read ``previous = rhs(1)`` to seed the M-M Carlson closure) +
`orpheus.numerics.iteration.KrylovAcceleration` (pre-fix
``preconditioner=None`` silently defaulted to invoking ``L.solve``
on each GMRES residual).

**Mechanism:** the curvilinear M-M closure
(Hébert §3.9.4 Eqs. 3.432-3.435) needs a half-angle ψ seed to start
the inward (``μ=−1``) sweep.  Pre-Phase-1.2 the seed came from the
lag-1 frame of the iterate via ``rhs(1)`` on the ``AngularFlux``
history machinery — adequate for the source-iteration loop where
each call passes a true ψ iterate with valid history.  GMRES feeds
**residual vectors** into the preconditioner; residuals have no
history (``rhs.history_depth == 0``); ``rhs(1)`` then returned the
``AngularFlux`` default (zeros via the cold-frame fallback).  The
M-M closure with a zero half-angle seed is a **poor preconditioner
for curvilinear** geometries — the Krylov polynomial
``M⁻¹·A·δ`` is non-identity-like and destabilises GMRES.  Slab
sweeps use no M-M closure (no curvilinear pole), so the slab
default-sweep precond worked correctly.  Sphere/cylinder
default-sweep precond exhibited 470× iteration-count blowup and
keff oscillation around the wrong fixed point.

**How it hid:**

- **Slab default-sweep worked.**  The M-M closure is curvilinear-
  only; slab sweeps never read the half-angle seed.  Step D's
  probe A (slab) converged in 1-2 outer iterations under either
  precond choice.  The diagnostic only surfaced when probe B
  exercised sphere.

- **The capability advertisement (``CAP_SOLVE``) was uniform across
  L+C operands**; the type system gave no signal that "this
  ``L.solve`` reads ``rhs(1)``" was a precondition any caller had
  to satisfy.  Stateful coupling was implicit in the source-
  iteration loop's design but not enforced anywhere.

- **The fallback path was implicit.**
  ``KrylovAcceleration(preconditioner=None)`` invoked
  ``self.L.solve(q)`` for every GMRES residual — there was no
  log, no warning, no type signal that the precondition (``q``
  must be a valid iterate with history) was being violated.
  Residual vectors had ``history_depth=0`` and ``rhs(1)`` silently
  returned the cold-frame default (zeros).

- **The slab cross-check ``passed`` for the wrong reason.**  The
  slab Krylov leg of every ``test_kinf_homogeneous`` case
  converged because slab doesn't trigger the M-M closure.  The
  sphere/cylinder Krylov legs failed; that was the diagnostic
  trigger.

**Empirical pre-fix evidence (R-1 Step D Probe B, 2026-05-19):**

```
sphere-2eg-krylov, 5 outer iterations:
  default (None) precond:  keff = 0.7..1.5..1.3..1.7..1.4  (oscillating)
  identity precond:        keff = 1.8750 (converges)

slab-2eg-krylov, 3 outer iterations:
  default (None) precond:  keff = 1.8750 (converges in 1 outer)
  identity precond:        keff = 1.8750 (converges in 1 outer)

Inner iteration count on sphere default-precond:  ~470× the count
on identity-precond (GMRES restarted repeatedly without progress).
```

**Empirical post-fix evidence (Phase 1.2, 2026-05-22):**

```
Phase 1.2 made StreamingCollisionOperator.solve a pure function:
  signature pre-fix:   solve(rhs)  →  read rhs(1) for Carlson seed
  signature post-fix:  solve(rhs, *, initial_guess=None)  →  seed
                       reads initial_guess explicitly; rhs(1) no
                       longer touched.

KrylovAcceleration(preconditioner=None) on a residual vector now
invokes LC.solve(residual) — initial_guess=None → cold-start M-M
seed (deterministic zero).  The M-M closure is still a poor
preconditioner for curvilinear under cold-start (numerical issue
tracked by #200), but the path is no longer silent and no longer
stateful.  Identity-precond is the production choice until #200
ships the block-inverse face preconditioner.

L1 kinf_homogeneous suite:           35 passed + 5 xfailed (unchanged)
L1 bridge-regression (Phase 1.3):   10/10 PASS
L1 precond-safety (Phase 1.3):       4/4 PASS
```

**Which test catches it:**

`tests/sn/test_krylov_curvilinear_precond_safety.py` (R-1 Step 4
Phase 1.3 promotion from
``derivations/diagnostics/diag_r1_step_d_probe_b_identity_precond.py``).
Two test functions, ``@pytest.mark.l1
@pytest.mark.catches("ERR-050") @pytest.mark.verifies(
"transport-cartesian", "sn-curvilinear-homogeneous-kinf-recovery")``:

1. ``test_identity_preconditioner_recovers_kinf`` — parametrised
   over {slab, sphere, cylinder}.  Pins the production contract
   that the typed Krylov inner solver with explicit identity
   preconditioner converges to analytical ``k_inf`` at rtol < 1e-8.

2. ``test_default_sweep_preconditioner_recovers_kinf_on_slab`` —
   slab only.  Pins the **structural-fix sentinel** for the
   silent-fallback bug class: slab default-sweep precond converges
   to ``k_inf`` because slab is the geometry where preconditioner
   quality (curvilinear M-M cold-start convergence) is NOT the
   confound.  If the structural fix is reverted (e.g. a stateful
   ``rhs(1)`` read re-introduced inside ``L.solve``), slab
   default-sweep would also destabilise on residual inputs because
   the residual carries no valid history.

Companion structural pin in
``tests/sn/test_invertible_operator.py::TestSolve::test_solve_forwards_explicit_initial_guess_to_sweep``
(``@pytest.mark.foundation``): spies on ``transport_sweep`` and
verifies the seed comes from the explicit ``initial_guess`` kwarg,
NOT from any ``rhs(1)`` history lookup.

**Why sphere/cylinder default-sweep is NOT pinned** by the L1
suite: the curvilinear M-M closure with cold-start zero seed is a
poor preconditioner for GMRES — a **numerical** issue, not the
structural bug-class fix.  Issue #200 designs the block-inverse
face preconditioner that restores sweep-as-preconditioner quality
on curvilinear; pinning a failing curvilinear default-sweep case
would lock in the wrong production state today.

**Lesson** (cf. `.claude/lessons.md` L19 + L21): default values for
behavioural parameters MUST either advertise their preconditions in
the type system (``CAP_STATELESS_INVERSE``-style capability gating)
OR require explicit caller choice (no default at all).  When you
catch a "default fell back silently to a stateful primitive" bug,
the **first** fix to consider is **structural unification**: are
the two consumers actually different applications of the same
operator?  If yes (here: sweep and matvec both apply the same
``(L+C)`` operator, only differing in which ψ they feed the M-M
closure as the seed), unify them through one strategy.  Reduce
strategies; don't add alternatives.  The structural fix (Phase
1.2) eliminated the bug class without needing a capability flag —
a stronger result than the original A2 plan because
``StreamingCollisionOperator.solve`` is now a pure function with no
silent paths to deprecate.

→ ERR-049 (convention drift) is the sibling defect — both arose
during the R-1 Step 4 carve when the operator-algebra layer met
the legacy sweep implementation; both closed by Phase 1.1 + Phase
1.2 of the consolidation plan.


---

## ERR-051 — `GalerkinProjection.assert_galerkin_idempotency` asserted :math:`\Pi R = I` instead of :math:`\Pi R = 4\pi \cdot I` (no-prefactor SH convention)

**Status:** **CAUGHT IN PHASE 1 P1.6 (2026-05-26, moment-space + layering plan).** Method deleted; sole caller deleted alongside.

**Failure mode:** **#6 (convention drift)** — the method's docstring named the discipline ("Galerkin idempotency :math:`\Pi R c = c`") but the no-prefactor SH basis convention installs a :math:`4\pi` factor that the addition-theorem reconstruction carries (:math:`R = (2\ell+1) \cdot S_0`, and the discrete Gram diagonal is :math:`4\pi/(2\ell+1)`, so the composition :math:`\Pi R = 4\pi \cdot I`, not :math:`I`). The discipline-statement was generic; the SH convention demanded a :math:`4\pi` correction the method did not apply.

**Date discovered:** 2026-05-26 Phase 1 audit while designing the V&V test suite for the ERR-039 endpoint (`tests/numerics/test_spherical_harmonic_space.py`).
**Module:** `orpheus.numerics.projection.GalerkinProjection.assert_galerkin_idempotency` (now deleted).

**Mechanism:** The method was added during Wave 0 alongside the `ProjectionOperator` ABC. The implementer drafted the assertion from the abstract Galerkin-idempotency statement ":math:`\Pi R = I` on the coefficient space", forgetting that the convention-bearing SH basis carries an implicit metric (the :math:`4\pi/(2\ell+1)` Gram diagonal) which modifies the discrete identity by a factor of :math:`4\pi`. The method called `np.testing.assert_allclose(self.apply(reconstruction.apply(c)), c)` — directly asserting :math:`\Pi R c = c`, not :math:`\Pi R c = 4\pi c`.

**How it hid:** **The method was never called against the production (Π, R) pair.** Its sole call site — `tests/numerics/test_projection_operators.py:368-381` (since deleted in P1.6) — deliberately built a **non-orthogonal Y matrix** so the method would raise:

```python
def test_method_signals_violation(self):
    # Construct a deliberately-broken Galerkin pair: use a
    # non-orthogonal Y matrix so Π R ≠ I.
    Y = rng.standard_normal((N, L + 1, 2 * L + 1))   # ← random Y, NOT a quadrature basis
    M = MomentProjection(weights=weights, Y=Y, L=L)
    R = HarmonicMomentReconstruction.from_Y(Y)
    with pytest.raises(AssertionError, match="Galerkin idempotency"):
        M.assert_galerkin_idempotency(R, c, atol=1e-10)
```

The test pattern is: "construct a broken pair → call the method → assert it raises." The method's docstring claim — that it would NOT raise on a correct pair — was never tested. Net effect: the method shipped with a wrong-identity error AND the test silently agreed because no one verified what "correct" should look like.

**How caught:** Phase 1 P1.5 designed `test_pi_R_is_4pi_identity_on_band_limited` for the ERR-039 endpoint. The test built the genuine :math:`(\Pi, R)` pair on a Lebedev-13 / :math:`L=3` quadrature and discovered :math:`\Pi R = 4\pi \cdot I`, NOT :math:`I`. Audit traced back to `assert_galerkin_idempotency`'s docstring; the discrepancy was the smoking gun. The sole caller's deliberately-broken-Y construction confirmed the test pattern was negative-only.

**Fix:** Delete the method and its sole caller (P1.6). The genuine identity is now pinned by:
- `tests/numerics/test_spherical_harmonic_space.py::test_pi_R_is_4pi_identity_on_band_limited` (uses `@pytest.mark.verifies("pi-r-equals-4pi-i")` — the existing Sphinx equation label).
- `tests/numerics/test_projection_operators.py::TestGalerkinIdempotencyOnLebedev::test_pi_R_is_identity_on_band_limited` (the sibling test that already pinned :math:`\Pi R = 4\pi \cdot I` correctly — pre-existed; was orphaned by the wrong-identity method's coexistence). *(2026-07-21 sync: this sibling file was retired by the Frame campaign; the surviving pin is the `test_spherical_harmonic_space.py` test above.)*

**Lesson:** **Every contract-validation method (``assert_X``, ``check_X``, ``verify_X``) MUST be tested against AT LEAST one correct instance where it must NOT raise AND AT LEAST one broken instance where it MUST raise.** Negative-only testing ("we fed it a broken case and it raised") validates the method's *raising behavior* but NOT its *invariant claim*. The bug surface in this defect was the TEST, not the method per se — the test never asked "does the method correctly distinguish correct from broken?", only "does it raise on broken?". Both halves are needed; either alone is a discipline failure that hides wrong-invariant claims indefinitely.

This generalizes the L11 lesson ("cross-checks must be structurally independent"): the cross-check here was not even procedural — it was *self-referential* (the broken Y was constructed precisely to make the wrong assertion succeed at raising). The structural-independence requirement applies to ALL test design, not just numerical cross-checks.

**Test reference:** `tests/numerics/test_spherical_harmonic_space.py::test_pi_R_is_4pi_identity_on_band_limited`. Tagged `@pytest.mark.l1` and `@pytest.mark.catches("ERR-051")` *(if the catches marker is added in a follow-up; currently the test is `@pytest.mark.catches("ERR-039")` because it was authored for the ERR-039 endpoint and ERR-051 surfaced as a side effect of writing it)*.

→ ERR-039 is the sibling: both are convention-drift failures of the SH Galerkin discipline that conflated the discrete metric with the abstract identity. ERR-039 was the operator-side conflation (`Π* = R`); ERR-051 was the verification-method-side conflation (`Π R = I`). Both closed by Phase 1 of the moment-space + layering plan (refactor/moment-space-and-layering, commits 0eb9cf3..c5be4b0).


## ERR-052 — Power iteration without per-step flux renormalisation: subcritical cases underflow to denormalised FP and converge to a meaningless keff ratio

**Status:** **CAUGHT 2026-05-26** while resolving pre-existing failures inherited from `refactor/sn-operator-algebra@62994ad` before starting Phase 3 of the moment-space + layering plan. Fixed in the same session.

**Failure mode:** **#3 (missing factor)** — the textbook power-iteration formula `ψ_{n+1} = (1/k_n) A^{-1} F ψ_n` requires an explicit renormalisation step `ψ_{n+1} /= ||ψ_{n+1}||` to prevent the iterate from growing (supercritical, `k > 1`) or decaying (subcritical, `k < 1`) geometrically. The codebase's two implementations — `orpheus.numerics.eigenvalue.power_iteration` (legacy) and `orpheus.numerics.iteration.KEigenvalue.solve` (canonical going forward, per P3.4 of the moment-space + layering plan) — **both** omitted the renormalisation.

**Date discovered:** 2026-05-26 during pre-Phase-3 baseline verification. The failing test (`tests/sn/test_boundary_conditions.py::TestSNBCSweepBehavior::test_vacuum_keff_lower_than_reflective`) was inherited unchanged from the `refactor/sn-operator-algebra` base; verified to fail identically there (commit 62994ad), confirming the bug pre-dates Phase 1.

**Module:** `orpheus.numerics.eigenvalue.power_iteration` (legacy, currently the `solve_sn` outer path) and `orpheus.numerics.iteration.KEigenvalue.solve` (canonical, will replace `power_iteration` in P3.4). Both have the same structural omission.

**Mechanism.** Inside the outer loop both implementations did:

```python
for n in range(1, max_iter + 1):
    flux_old = flux.copy()
    keff_old = keff
    q = solver.compute_fission_source(flux, keff)   # F·ψ / k
    flux = solver.solve_fixed_source(q, flux)       # (L+C-S)^{-1} · q
    keff = solver.compute_keff(flux)                # ratio of integrals
    if solver.converged(...): break
```

Without renormalisation, the iterate scales by approximately `keff` per step (the power-method ratio). For a 2cm slab with the 2-group SN benchmark material (`k_inf = 1.875`):

- **Reflective BCs** (no leakage, physically `k = k_inf = 1.875 > 1`): flux grows ~1.875× per outer step, but converges in 3-4 outers — the growth never reaches FP overflow. Test passes **accidentally**.
- **Vacuum BCs** (the slab is supercritical with `k ≈ 1.67 < k_inf`, but the test only needs `k_vac < k_refl`, which is the physically correct ordering): flux decays each step if the leakage is enough to push the operator's spectral radius below unity in the iterate's working magnitude. After ~30-60 outer steps the flux saturates at IEEE-754 **denormalised** magnitude (`~1e-39` to `1e-43`); both `compute_keff = ∫νΣ_f φ / ∫Σ_a φ` and the convergence-test residual become `0/0` ratios that round to spurious values. The legacy `converged()` method's `max(||ψ||, 1e-30)` floor magnifies the false-positive (the floor caps the denominator and produces a deceptively-small `dphi`); the iteration claims `converged=True` at a meaningless fixed point. The two inner solvers (`source_iteration` vs `krylov`) project onto different residual modes during the collapse, hence each reports a different bogus keff (SI: 1.94 ≈ k_g1_pure = 2.00 — degenerated to the thermal-only spectral shape; Krylov: 1.67 — close to the physical answer but with denormalised flux).

The keff ratio `(F·φ, φ) / (A·φ, φ)` is scale-invariant in φ, so per-step renormalisation `ψ /= ||ψ||` preserves keff while keeping the iterate at unit norm — the textbook power-iteration form. Post-fix all paths converge to **keff = 1.6693** (slab is subcritically supercritical: k_inf=1.875 reduced by leakage to 1.67, still > 1), with **psi_max = 0.089** (well-conditioned, not denormalised), in **n_outer = 6** (not the 64-step cap), with **SI and Krylov agreeing** to 1e-9.

**How it hid.** Every existing L1 eigenvalue test in `tests/sn/l1_analytical/test_kinf_homogeneous.py` (and the in-suite `test_invertible_operator.py` regression matrix) uses **reflective** BCs. Reflective always gives `keff = k_inf > 1`, where the flux grows but the iteration terminates before the growth blows up. The vacuum-eigenvalue path had **zero L1 coverage**; the only test exercising vacuum + power iteration was `test_vacuum_keff_lower_than_reflective`, which has been failing on the base branch since at least the `refactor/sn-operator-algebra` cut (verified at commit 62994ad). The failure pre-existed Phase 1 — Phase 1 inherited and did not introduce it.

Compounding the gap: `solve_sn` (`orpheus/sn/solver.py:992-996`) hardcodes `IterationHistory(..., converged=True)` regardless of the solver's actual `converged()` return value, masking the convergence-state defect from any caller that only inspects the history.

**How caught.** Pre-Phase-3 baseline run: progressive-chunk pytest exposed 2 failures on the moment-space-and-layering branch. The user (per `feedback_fix_bugs_immediately`) directed fix-on-this-branch-before-Phase-3. Static analysis by `numerics-investigator` narrowed the suspect set; probe-cascade discrimination (`source_iteration` vs `krylov`, isolated subprocesses, `python -c` vs `python script.py` path differences) pinned the failure to the iteration loop; reading `power_iteration` source revealed the missing `flux /= ||flux||`.

**Fix.** Add the standard power-method renormalisation at the end of each outer step, before the keff update:

```python
norm = float(np.linalg.norm(flux))
if norm > 0.0:
    flux = flux / norm
keff = solver.compute_keff(flux)
```

Applied identically in:
- `orpheus/numerics/eigenvalue.py:power_iteration` (legacy; retires in P3.4)
- `orpheus/numerics/iteration.py:KEigenvalue.solve` (canonical going forward)

Both call sites annotated with an `ERR-052` reference in the inline comment so future readers see the load-bearing reason.

**Lesson.** **Every BC type needs its own L1 eigenvalue test at multi-group.** "Reflective passes" does NOT generalise — reflective is the easy case for power iteration (flux growth is bounded by the small iteration count needed for `k > 1` cases); vacuum, white, and albedo BCs exercise the subcritical regime where the growth ratio is `< 1` and the iterate decays. Coverage by BC type is a separate axis from coverage by mesh / multigroup / scattering order. The moment-space + layering plan's P3.4 verification programme will install an `L1` test matrix indexed by `(BC, geometry, group count)` — that matrix would have caught ERR-052 the first time vacuum BCs were exercised on a multiplying medium.

**Secondary surface — Krylov inner-iteration budget.** The unit-production-rate normalisation alters the GMRES initial guess at each outer step (the previous un-normalised trajectory inherited a warmed-up subspace; the normalised trajectory does not). On `sphere-2eg-krylov` this raised the per-outer-iter GMRES count from ~50 to ~600 — and the L1 test `_TIGHT_KW` budget of `max_inner=300` was no longer sufficient. The inner solve was hitting the cap, returning an under-converged result, and the outer iteration accumulated ~2.4e-7 keff drift before claiming convergence. The fix was `max_inner=300 → 1000` in `tests/sn/l1_analytical/test_kinf_homogeneous.py::_TIGHT_KW`, restoring FP-precision keff for all 28 `(coord, ng, inner_solver)` variants. Issue #200 (block-inverse preconditioner for Krylov on the typed AngularFlux algebra) tracks the longer-term reduction. The production `solve_sn` default (`max_inner=200`) is unaffected — this is purely an L1 verification budget for the tightest reference-comparison gate.

Secondary lesson: **`IterationHistory.converged=True` was hardcoded** in `orpheus/sn/solver.py:992-996`, decoupled from the solver's actual convergence flag. This is a latent bug worth a follow-up fix (low severity — the keff value is correct post-ERR-052; only the `converged` field is misleading). Tracked as P3.4 close-out work.

**Test reference:** `tests/sn/test_boundary_conditions.py::TestSNBCSweepBehavior::test_vacuum_keff_lower_than_reflective` — tagged `@pytest.mark.catches("ERR-052")`. Companion diagnostic: `derivations/diagnostics/diag_vacuum_bc_eigenvalue_divergence.py` runs the three discriminating probes (reflective baseline, vacuum + SI, vacuum + Krylov) and confirms SI/Krylov agreement post-fix. The diag script also documents a session-level gotcha — standalone scripts under `derivations/diagnostics/` must prepend the repo root to `sys.path` to load the worktree's `orpheus`; otherwise the venv's `pip install -e .` silently resolves to the main checkout, giving stale-fix false-negatives.

→ Probe path: see `.claude/agent-memory/numerics-investigator/vacuum_bc_eigenvalue_divergence.md` for the full hypothesis cascade.

## ERR-053 — Hardcoded GMRES `restart=min(50, full_size)` clamp + discarded scipy `info` flag silently truncate the Krylov subspace and consume the unconverged iterate as the inverse

**Status:** **CAUGHT 2026-05-28** during the D-H.1 (TimedFullField composite carrier) trunk migration close-out — surfaced by `tests/sn/spatial/test_sweep_vs_apply_consistency.py::test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere` failing in the post-Stage-4 envelope. The bug pre-dates D-H.1 entirely; verified to fail identically on commit `9d02ade` (the D-G consolidation, well before D-H.1b.1).

**Failure mode:** **#3 (missing factor)** + **#7 (test-design failure / MMS simplification bias)** in superposition. The structural defect is a SUBSPACE-DIMENSION undershoot at the GMRES call site; the test-coverage defect is that no existing L1 anchor probed a homogeneous-reflective Krylov-eigenvalue case at `n_unknowns > 50` with a strict `keff_tol`. The two together let the bug ship under "all gates green" while quietly producing keff errors up to **47.3%** on curvilinear meshes.

**Date discovered:** 2026-05-28. The failure surfaced when the D-H.1c stage-4 envelope expanded to include the `tests/sn/spatial/` suite (previously excluded for wall-clock per L16). The numerics-investigator promotion + restart-sweep diagnostic pinned the root cause.

**Module:** `orpheus.sn.solver` (lines 716 and 1477 — the `_solve_krylov` method on `SNSolver` and the module-level `_solve_curvilinear_krylov` helper) and `orpheus.numerics.iteration.KrylovAcceleration.solve` (lines 735, 746 — both scipy.gmres call branches).

**Mechanism.** Two defects in superposition, mutually concealing.

**Defect 1 — Subspace truncation at the call site** (`orpheus/sn/solver.py:716`, `:1477`):

```python
krylov = KrylovAcceleration(
    LC, scattering_op, ZeroOperator(),
    preconditioner=lambda q: q,
    tol=inner_tol, max_iter=max_inner,
    restart=min(50, N * ng * nx * ny),   # ← THE BUG
)
```

GMRES with `restart=R` truncates the Krylov subspace to dimension `R`; any direction outside that subspace cannot be represented in the iterate, regardless of `maxiter`. The clamp `min(50, full_size)` is dimensionless intent (presumably "don't waste memory on tiny problems") but the threshold `50` is small enough that EVERY non-trivial curvilinear problem trips it: a sphere with N=8 ordinates × 2 groups × 20 cells × 1 ny = 320 unknowns, plus 8 boundary trace ordinates × 2 groups = 16 face slots = **328 unknowns total**. Restart=50 vs full=328 ⇒ GMRES is structurally unable to reach the true `(L − S − F)⁻¹ q` from any non-trivial initial guess.

**Defect 2 — Discarded convergence flag** (`orpheus/numerics/iteration.py:735` and `:746`, both scipy branches):

```python
solution, _info = spla.gmres(...)
```

scipy.gmres returns `(solution, info)` where `info > 0` signals "did not converge within `maxiter`" and the returned `solution` is the best-effort iterate. The underscore-prefix discards the flag. Without the warning, `KrylovAcceleration.solve` returns the unconverged iterate as if it were the inverse; the consumer (the eigenvalue outer loop) computes a `keff` from a flux that doesn't satisfy `(L − S − F) ψ = q_ext` even to its claimed tolerance.

**Why "tighter inner_tol" did NOT rescue the test** (the misleading first hypothesis): scipy's `rtol` controls the EXIT criterion (relative residual < rtol), not the SUBSPACE dimension. With restart=50, no value of `rtol` can produce convergence — the subspace simply doesn't contain the solution direction. The investigator's repro (`/tmp/repro_krylov_sphere.py`) showed keff = 1.4019239124 identically for `inner_tol ∈ {1e-8, 1e-10, 1e-12}` — a clean tolerance-isn't-it falsification.

**Error growth with refinement** (the load-bearing diagnostic): mesh-refinement table from the investigator's step-5 diagnostic (`derivations/diagnostics/diag_krylov_si_homogeneous_sphere_step5_mesh_scaling.py`):

```
n_cells   SI_keff         SI_err     KR_keff         KR_err
   5      1.8750000000   1.069e-11   1.8750000004   4.103e-10
   10     1.8750000000   1.097e-11   1.9152507886   4.025e-02
   16     1.8750000000   1.104e-11   1.5987097760   2.763e-01
   20     1.8750000000   1.106e-11   1.4019239124   4.731e-01
```

The signature is unmistakable: SI is correct at machine precision across the refinement series; Krylov is correct only when `n_unknowns ≤ 50` (the n_cells=5 row), then diverges as the mesh refines and the natural subspace dimension grows past the hardcoded clamp. **Error growing with refinement** is the canonical structural-defect signature — a tolerance issue would produce uniform error or monotone-decreasing error with refinement, not divergence.

**Restart sweep — direct scipy** (step-8a):

```
restart=  50  info=200  true_||r||=1.273e+00         ← unconverged, info ignored upstream
restart= 100  info=  0  true_||r||=1.421e-12
restart= 200  info=  0  true_||r||=1.421e-12
restart= 320  info=  0  true_||r||=1.251e-12
```

`info=200` at restart=50 == "did not converge in 200 outer iterations of the restarted scheme"; the iterate's true residual is **O(1)**, not O(rtol). The discarded `info` flag (Defect 2) is what let this slip into the keff calculation.

**End-to-end fix verification** (step-8b — manual outer loop, no production-rate rescale, no `_psi_typed` cache):

```
restart  keff           err         last_res
   50    1.5498783237   3.251e-01    5.677e+00   ← production hardcode
  100    1.8750000000   6.395e-14    9.193e-13   ← bug dissolves
  200    1.8750000000   6.972e-14    9.988e-13
```

Once the subspace is large enough to span the solution direction, GMRES converges to floating-point precision.

**How it hid.** Three reasons:

1. **L1 anchor uses tight `inner_tol=1e-12` AND a small mesh** (`tests/sn/test_krylov_curvilinear_precond_safety.py::test_identity_preconditioner_recovers_kinf`, n_cells=10). For n_cells=10, the natural subspace happens to fit within the 50-dimension clamp on the dominant eigenmode for this test's specific operator; the iteration coincidentally projects correctly. Mesh refinement past n_cells≈16 was not in the L1 matrix.

2. **The failing test (`test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere`) was a `tests/sn/spatial/` member**, and the spatial suite was excluded from the leaf envelope during D-H.1b/c (per L16, wall-clock constraint — full-suite execution was hitting 70 minutes). The bug had been failing on every base commit since `refactor/sn-operator-algebra@62994ad` was cut, but no agent had run the spatial suite in any prior session.

3. **Discarded `info` flag silently consumes the unconverged iterate**. Without Defect 2, scipy would have surfaced the failure as a `LinAlgWarning` at the consumer; with Defect 2, the failure is invisible to every consumer including the eigenvalue outer loop.

The interaction of #1+#2+#3 is the test-design failure (Mode #7 — the existing L1 coverage's parameter range did not stress-test the bug-rich corner of the parameter space).

**How caught.** Post-Stage-4 of D-H.1c (the trunk migration to TimedFullField composite carrier), the leaf-envelope retest expanded to include `tests/sn/spatial/` for the first time. `test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere` failed with `SI keff = 1.875, Krylov keff = 1.4019` — a 25% drift that was immediately suspicious (vs. the L1 anchor's rtol < 1e-8 result on the same operator). Pre-D-H.1 bisection on commit `9d02ade` confirmed the bug pre-dated the migration. Static analysis surfaced the `restart=min(50, ...)` clamp; numerics-investigator's step-8 restart-sweep diagnostic pinned the mechanism by direct scipy invocation, isolating the subspace truncation from the `info` discard.

**Fix.** Two single-line changes preserving the public-API contract (`inner_tol`, `max_inner`, `keff_tol`, `flux_tol` semantics unchanged):

1. `orpheus/sn/solver.py:716` (`SNSolver._solve_krylov`) and `:1477` (`_solve_curvilinear_krylov` module helper):

```python
# Before:
restart=min(50, N * ng * nx * ny),
# After:
restart=N * ng * nx * ny,
```

Restart at the full problem size disables subspace truncation; GMRES converges from any well-defined initial guess.

2. `orpheus/numerics/iteration.py:735, 746` (both scipy branches in `KrylovAcceleration.solve`):

```python
# Before:
solution, _info = spla.gmres(...)
# After:
solution, info = spla.gmres(...)
if info != 0:
    warnings.warn(
        f"KrylovAcceleration.solve: scipy.sparse.linalg.gmres returned info={info} "
        f"(not converged within maxiter={self.max_iter}; restart={...}; rtol={self.tol}).  "
        f"Returning best-effort iterate; residual_history tail = ...  "
        f"Tighten ``restart`` to ``n`` (full size) if the Krylov subspace is being "
        f"truncated; see ERR-053.",
        RuntimeWarning, stacklevel=2,
    )
```

Surfaces non-convergence so the consumer can choose to retry / raise / proceed. RuntimeWarning rather than raise preserves call-site compatibility for callers that legitimately tolerate non-convergence (e.g., diagnostic probes).

**Verification.** Post-fix:

```
inner_tol=1e-08: keff=1.8749999997     ← was 1.4019239124
inner_tol=1e-10: keff=1.8750000000
inner_tol=1e-12: keff=1.8750000000
```

`test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere` — 1 passed in 2.32s (was failing with `abs(diff) = 4.73e-01 > 1e-6`).

**Lesson.** **Subspace-dimension caps are SILENT failure modes for iterative linear solvers.** Unlike tolerance caps (which fail with a residual signal), subspace truncation produces a structurally-wrong answer with no observable signal at the call site. The compounding bug (`_info` discard) created a silent failure of a silent failure. The defense is twofold: (a) NEVER discard the convergence flag of a scipy iterative solver — promote it to at least a warning; (b) NEVER hardcode a subspace size below the natural problem dimension without an explicit MAX_INNER-shape verification gate.

**Secondary defense at the test level**: mesh-refinement convergence is the canonical structural signature for distinguishing tolerance defects from subspace-dimension defects. Tolerance defects produce uniform or monotone-decreasing error with refinement; subspace defects produce DIVERGING error with refinement (because the natural subspace grows past the cap). The new permanent regression test (`tests/sn/test_krylov_restart_signature.py`, promoted from `diag_krylov_si_homogeneous_sphere_step5_mesh_scaling.py` per the investigator's recommendation) pins this signature at L1.

**Test reference:** `tests/sn/spatial/test_sweep_vs_apply_consistency.py::test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere` (existing — pinned by inheritance), plus the new mesh-refinement regression catcher under `tests/sn/` (this commit), plus the restart-sweep direct-scipy diagnostic that confirms `info` discard at the kernel boundary. All three carry `@pytest.mark.catches("ERR-053")`.

→ Probe path: see the 8 diagnostic scripts at `derivations/diagnostics/diag_krylov_si_homogeneous_sphere_step{1..8}_*.py` for the full bisection cascade.


## ERR-054 — `ordinate_scan` Blelloch closed-form `cumprod_a · (psi_0 + cumsum(b/cumprod_a))` produces NaN when any chain entry of `a` is exactly 0 (cylindrical pole-cell algebraic resonance)

**Date:** 2026-05-29
**Module:** `sn` (`orpheus/sn/spatial/scan.py`)
**GitHub issue:** [#209](https://github.com/deOliveira-R/ORPHEUS/issues/209)
**Failure mode:** #4 (wrong recursion — closed-form numerically catastrophic in a documented but unenforced regime)

**Symptom.** `solve_sn(inner_solver="source_iteration")` on a 1-D cylindrical reflective mesh at the EXACT configuration `(thickness=2.0 cm, n_cells=20, LevelSymmetric S8, mixture A 2G)` returns `keff = NaN` and times out at 475 s under default `max_outer=500 × max_inner=500` caps. The same problem solved via `inner_solver="krylov"` converges to `keff = 1.875` (the analytical `k_inf = νΣ_f / Σ_a`) in 0.03 s. Neither `n_cells = 19` nor `n_cells = 21` triggers the failure; only the sharp algebraic resonance at the exact `(thick=2, n=20, μ_x=-1/√20, Σ_t=1)` point.

**Root cause.** `orpheus/sn/spatial/scan.py:138` evaluates the per-ordinate first-order linear recurrence `ψ[i+1] = a[i]·ψ[i] + b[i]` via the Blelloch closed form

```python
cumprod_a = np.cumprod(a, axis=0)
return cumprod_a * (psi_0 + np.cumsum(b / cumprod_a, axis=0))
```

This is correct in real arithmetic but numerically catastrophic in IEEE-754 when any chain entry `a[i] = 0` exactly. The cumprod then collapses to zero from cell `i` onward; the divide `b / cumprod_a` produces `±Inf` or NaN; `np.cumsum` propagates the NaN forward; and the final `cumprod_a · (...) = 0 · NaN = NaN`, not zero as the math would suggest. The mathematically equivalent explicit recurrence is FINITE at the same chain: at `a[i] = 0`, the per-cell update degenerates to `ψ[i+1] = b[i]`, which is a physically meaningful "fully attenuated" exit flux.

The cache (`orpheus/sn/spatial/sweep_cache.py:444-449`) produces `a` via

```python
a = 2|μ| · A_total / (2|μ| · A_down + dA_w · c_out + Σ_t · V) − 1
```

At the **cylindrical pole cell** the inner radial face has zero area, so `A_down = 0`. The formula reduces to `a = 2|μ|·A_total / (dA_w·c_out + Σ_t·V) − 1`, and `a = 0` ⇔ `2|μ|·A_total = dA_w·c_out + Σ_t·V`. At `μ_x = -1/√20` (smallest |μ| in LS-S8), `dr = 0.1`, `Σ_t = 1.0` (mixture A group 1), the identity holds **bit-exactly** at `2π/5 = 1.2566370614359172`. `cache.a_attenuation[ord=72, g=1, chain_pos=19] == 0.0` is the smoking-gun cache entry.

**How it hid.**

1. **Existing scan-form tests cover the wrong regime.** `tests/sn/spatial/test_ordinate_scan.py::test_ordinate_scan_small_attenuation` uses `a ∈ [0.05, 0.2]` (positive, bounded away from 0). `test_ordinate_scan_zero_attenuation` uses `a ≡ 1` (the opposite limit). No existing test covers `a ∋ 0` — the exact pole-cell pathology.

2. **The docstring acknowledged the regime but did not enforce it.** `scan.py` lines 126–135 explicitly noted "requires `cumprod_a` to stay finite and bounded away from zero ... For `a → 0` ... outside DD's normal operating envelope, consult the test catalog". The caveat lived in prose; no positive contract test pinned it. Anti-pattern #10 in `vv-principles` ("docstring caveat without enforcement").

3. **Production SI cylindrical tests miss the resonance.** Integration suites use `n_cells ∈ {40, 80, 160}` (dr ∈ {0.05, 0.025, 0.0125}); none lands on `dr = 0.1` with `μ_x = -1/√20, Σ_t = 1`. The L1 standoff suite (cited in the bug report as "slowed by the bug") uses `thickness=2.0, n_cells=40` and does NOT actually hit this resonance.

4. **1-G eigenvalue degeneracy + homogeneous-reflective shape invariance.** Even when the bug is active, `k = νΣ_f / Σ_a` is independent of the angular flux shape; a converged SI returns `k_inf` regardless of internal redistribution errors. NaN is the ONE failure that cannot be smuggled past the eigenvalue identity — it propagates.

5. **Krylov bypasses the buggy code path entirely.** `transport_operator_matvec_unified` and the per-geometry matvec helpers in `orpheus/sn/operator.py` do NOT import `ordinate_scan`; only `_sweep_1d_unified` (the SI sweep path) does. The structural-independence assertion is verified in `diag_si_cyl_20cell_nan_step5_root_cause.py::test_krylov_avoids_ordinate_scan_path`.

**Which test catches it.** Permanent regression catcher: `tests/sn/test_si_cyl_20cell_nan_regression.py`. Pre-fix this test FAILS on:

* `test_si_returns_finite_keff` (SI returns NaN — the bug class signature).
* `test_ordinate_scan_at_a_zero_returns_finite_via_loop` (the scan-form contract test, structurally independent of any solver).

Post-fix (when the Blelloch closed form is replaced with a numerically-stable pair-monoid scan), both pass; `test_si_agrees_with_kinf_at_resonance` (currently `@pytest.mark.slow`) becomes the L1 correctness pin.

Diagnostic scripts:

* `derivations/diagnostics/diag_si_cyl_20cell_nan_step1_characterize.py` — 6 tests pinning the sharp-resonance fingerprint (`n_cells = 20` only; Krylov works on the same problem).
* `derivations/diagnostics/diag_si_cyl_20cell_nan_step5_root_cause.py` — 4 tests pinning the cache-level `a = 0` algebraic identity, the `ordinate_scan` NaN at the failing chain, the explicit-loop finiteness, and the Krylov-bypass structural invariant.

When the fix lands, both regression tests carry `@pytest.mark.catches("ERR-054")`.

**Fix family.** Replace the Blelloch closed form with a numerically-stable Blelloch variant. Three viable options:

1. **Pair-monoid prefix scan.** Compose `(α, β) ⊕ (α', β') = (α·α', α'·β + β')` via an explicit associative prefix scan. The existing `tests/sn/spatial/test_ordinate_scan.py::test_pair_monoid_associativity` already verifies the algebra. No division anywhere; vectorises across `(K, ng)` identically. **Preferred — cleanest path.**

2. **Fallback to explicit loop** at chain cells where `|a| < ε`. Hybrid; loses the all-numpy uniformity.

3. **Brent blocked scan** with bounded condition number per block. Over-engineered for 1-D SN.

Coordination: the fix lives in `orpheus/sn/spatial/scan.py:80-138`; it is orthogonal to the D-H.2-C5 `angular_flux.py` retirement. It can land independently or alongside.

**Lesson.**

1. **Closed-form algorithms with regime caveats need positive contract tests at the regime boundary.** A docstring claiming "stable when `a` is bounded away from 0" is NOT a contract; only an executable test that constructs a chain with `a = 0` and asserts finiteness is. ERR-054 demonstrates that documented-but-untested regime limits are silent-failure waiting to happen.

2. **`vv-principles` Anti-pattern #10 ("convergence rate is correct" ≠ "result is correct") generalises to algorithmic stability:** a closed-form's correctness in real arithmetic is not the same as its correctness in IEEE-754. Real-arithmetic equivalence is not numerical equivalence. The pair-monoid SCAN proof on paper does not save the IEEE-754 cumprod_a from collapsing to zero.

3. **Sharp algebraic resonances (single-point failures) are an under-tested bug class.** Mesh-refinement convergence tests (n=10, 40, 80, 160) implicitly assume the bug is a smooth function of mesh; but a clean algebraic identity like `2|μ|·A_total = dA_w·c_out + Σ_t·V` is single-point and can be missed by any refinement sweep that doesn't probe the exact point. The defense is dimensional-analysis-derived corner probes: for every closed-form `f / g − 1` numerator-denominator structure, write a test that constructs the exact (g, f) coincidence and asserts the resulting value is finite. Catalogue this for follow-up: there is no `numerical-bug-signatures` Signature 8 for "closed-form-divides-by-zero-at-algebraic-resonance"; the catalog gap is itself a finding.

4. **Krylov-versus-SI structural divergence is a load-bearing cross-check.** When two solver paths share the same operator construction but differ in execution algorithm, agreement is information; disagreement at the SAME problem is a diagnostic localiser. Issue #209 was localised in <2 hours because Krylov-vs-SI disagreement was already in the user's empirical table.

**Test reference:** `tests/sn/test_si_cyl_20cell_nan_regression.py` (this commit), with `@pytest.mark.catches("ERR-054")` to be added when the fix lands.

→ Probe path: see `derivations/diagnostics/diag_si_cyl_20cell_nan_step{1,5}_*.py` for the cascade. The cascade is two-step (no need for step 2/3/4 isolation because the failing path was named directly by the FP-warning traceback at step 1); the methodology is a degenerate case of the standard 8-step cascade where step 1's traceback short-circuits the isolation.

## ERR-055 — Curvilinear sweep regression tests fed `sig_t` / `Q` in the obsolete `(nx, ng, ny)` layout after the production contract flipped to PR-INDEX-5 `(ng, nx, ny)` / `(N, ng, nx, ny)` — `IndexError` at `CollisionCache.from_geometry`

**Date:** 2026-06-01
**Module:** `sn` / `tests` (`tests/sn/test_spherical.py`, `tests/sn/test_cylindrical.py`; crash surfaces at `orpheus/sn/spatial/sweep_cache.py:431`)
**Failure mode:** #6 (convention drift — definition site vs usage site disagree), residing in the TEST fixtures, NOT in production code.

**Symptom.** Six SN curvilinear sweep-regression tests crash with `IndexError: index 9 is out of bounds for axis 1 with size 1` at `sweep_cache.py:431` (`sig_t[:, geom.chain_idx].transpose(1, 0, 2)`). Affected: `test_spherical.py::TestSphericalSweepRegression::{test_uniform_source_converges_to_Q_over_sigt, test_single_sweep_all_finite}`, `::TestMultiGroupMultiRegionSpherical::test_fixed_source_flux_bounded`, `test_cylindrical.py::TestCylindricalSweepRegression::test_single_sweep_all_finite`, `::TestMultiGroupMultiRegion::{test_redistribution_telescoping_conservation, test_single_cell_uniform_source_equilibrium}`.

**Root cause.** Each test calls the internal `_sweep_1d_unified(Q, sig_t, sn_mesh, boundary_flux)` directly with bare ndarrays in the **obsolete** `(nx, ng, ny)` layout — e.g. `sig_t = np.full((10, 1, 1), 0.5)`, `Q = np.ones((10, 1, 1))` for a 10-cell ng=1 mesh. The production contract (`sweep.py:167-168`, `sweep.py:164`) is `sig_t (ng, nx, ny)` and source `(N, ng, nx, ny)`. With ng=1 these are NOT interchangeable: `_ensure_coll_cache` slices `sig_t[:, :, 0]` → `(10, 1)`, and `from_geometry` reads axis 0 as `ng=10`, axis 1 as `nx=1`; `geom.chain_idx` then indexes the 10-cell chain into a size-1 axis → `IndexError`. The bug was introduced when commit `6cfdfd4` (Issue #196 PR-INDEX-2) flipped `CollisionCache` to the `(ng, nx)` layout and the principled `(N, ng, nx, ny)` source / `(ng, nx, ny)` sig_t contract landed; the production producers (`solver.py`, `operator.py`, `material_xs_field.py:372` `xs.sig_t.T.reshape(ng, nx, ny)`) were migrated but these six direct-call tests were not. `_sweep_1d_unified` is **production** (the sole 1-D sweep body, reached by `transport_sweep` from `solver.py:970` and `operator.py:1946`); it is NOT a dying path — so the principled fix is a test migration, not a helper retirement.

**How it hid.**

1. **The crash is deterministic but only on the DIRECT-call tests.** Every production sweep flows through `transport_sweep` → `PerOrdinateSource.from_isotropic(..., sn_mesh)`, which builds the source in the principled `(N, ng, nx, ny)` layout by construction. Only ad-hoc tests that bypass the producer and hand-build bare arrays could carry the obsolete layout — these six did. The green matvec twin (`test_unified_matvec_cylinder.py`) constructs `sigma_t = np.full((ng, n_cells, 1), 2.0)` (correct layout) and never crashed, masking the test-side drift.

2. **The degenerate ng=1 layout aliasing.** `(nx=10, ng=1, ny=1)` and `(ng=1, nx=10, ny=1)` have the same total element count and look superficially similar; only the trailing-slice + chain-index step exposes the axis swap. An ng≥2 test would have failed louder/earlier (axis-1 size 2 vs nx mismatch), but these six are all ng=1.

3. **Pattern 7 violation in the test layer.** The tests re-applied the array-layout convention at the consumer (hand-built bare arrays) instead of routing through the producer. A future layout flip re-opens the drift silently because nothing pins the test inputs to the production convention.

**Which test catches it.** The fix migrates all six tests onto the production producer: `sig_t` shaped `(ng, nx, ny)` and source built via `PerOrdinateSource.from_isotropic(Q_iso(ng,nx,ny), sn_mesh)`, swept through `transport_sweep`. Output indexing updated `phi[:, 0, 0] → phi[0, :, 0]` and angular `ang[n, :, 0, 0] → ang[n, 0, :, 0]` to match `(ng, nx, ny)` / `(N, ng, nx, ny)`. The six migrated tests ARE the permanent regression gates (they now exercise the same convention production obeys; a future re-flip fails them).

**Lesson.**

1. **`coding-elegance` Pattern 7 applies to test fixtures, not just production.** When a test bypasses the production producer to hand-build inputs, it re-applies a layout/normalisation convention at the consumer — an independent opportunity for drift. Route tests through the SAME producer the solver uses (`PerOrdinateSource.from_isotropic`, `transport_sweep`) so layout drift cannot silently re-open. A test that calls an internal `_helper` directly with bare arrays is a Pattern-7 landmine.

2. **A layout/convention migration is incomplete until the direct-call tests are migrated too.** PR-INDEX-2/PR-INDEX-5 migrated the production producers but left six direct `_sweep_1d_unified` callers on the old layout. Per `feedback_retirement_means_test_migration`, flipping a convention means rewiring its tests to the new convention — and the test inputs must match the producer the new code expects.

3. **ng=1 degeneracy hides axis-swap convention drift.** `(nx, ng, ny)` with ng=1 aliases against `(ng, nx, ny)` in element count and shape-rank; the swap only declares itself at an index/slice that distinguishes the two axes. ≥2 groups would have surfaced it sooner — another instance of cross-cutting hygiene rule H1 (1-group degeneracy) operating at the data-layout level.

**Test reference:** `tests/sn/test_spherical.py::TestSphericalSweepRegression::{test_uniform_source_converges_to_Q_over_sigt, test_single_sweep_all_finite}`, `::TestMultiGroupMultiRegionSpherical::test_fixed_source_flux_bounded`, `tests/sn/test_cylindrical.py::TestCylindricalSweepRegression::test_single_sweep_all_finite`, `::TestMultiGroupMultiRegion::{test_redistribution_telescoping_conservation, test_single_cell_uniform_source_equilibrium}` (all six migrated this commit). No `@pytest.mark.catches` needed — the bug was in the test fixtures, not a production code path; the migrated tests gate the production convention by routing through the producer.

---

## ERR-056 — Octant-group Gauss-Seidel schedule reflected a boundary face after only the FIRST octant group outflowing it, absorbing the not-yet-swept octants' SEED value instead of their outflow → wrong fixed point on diagonal/spherical cubatures

**Date:** 2026-06-05
**Module:** `sn` (`orpheus/sn/sweep_schedule.py::SweepSchedule.gauss_seidel`; surfaces via `orpheus/sn/solver.py::_GaussSeidelResolvent` on the 2-D fixed-source SI path)
**Failure mode:** a scheduling / data-dependency bug (premature read of stale ephemeral data) — adjacent to mode #4 (wrong recursion / index drift), but at the octant-SCHEDULE level rather than the cell recurrence.

**Symptom.** With the Phase 3 boundary Gauss-Seidel default (`inner_schedule="gauss_seidel"`), two 2-D fixed-source tests using a `Quadrature.lebedev(order=17)` cubature converged to the WRONG flux: `test_2d_octant_sweep_closed_form_anchor` (all-reflective uniform → should be flat `5.882`) gave a non-flat `~3.4` (max |diff| 13.6); `test_2d_heterogeneous_si_krylov_equivalence` (fuel|moderator, vacuum-x) gave a FLAT `2.748` where Krylov gave the correct spatially-varying `2.43/2.40/...` (max rel-diff 1.23). The G-S iteration CONVERGED (the iterate stabilised) — to a wrong fixed point.

**Root cause.** `SweepSchedule.gauss_seidel` originally assigned each in-plane octant group its OWN outgoing reflective faces (`reflect_faces = _outgoing_faces(label) ∩ reflective`), reflecting each face right after its octant group swept. On an AXIS-ALIGNED quadrature (`product` — single-face octants) each reflective face is outflowed by exactly one octant group, so this is correct. But on a DIAGONAL / spherical cubature (`lebedev`, `level_symmetric` — each octant outflows TWO faces) a face is SHARED by ≥2 octant groups (e.g. `xmax` ← every +x octant: `(+1,0)`, `(+1,+1)`, `(+1,−1)`). Reflecting `xmax` after only the FIRST +x group absorbed the wavefront's `xmax` edge while the OTHER +x octants were still unswept — and because the `WavefrontFlux` is rebuilt + ι_*-seeded each `.solve`, those unswept octants' `xmax`-outflow slots still held the INFLOW SEED, not real outflow. The reflect (`R · seed`) produced garbage inflow, and the iteration converged to a fixed point of the WRONG operator. The seed contamination does NOT self-correct at convergence (unlike a stale-but-real lag) because the ephemeral wavefront never contains the unswept octants' outflow until they are swept.

**How it hid.**

1. **Axis-aligned quadratures are immune.** The first-built diagnostics + the production verify all used `product(n_mu=2, n_phi=4)` (4 single-face in-plane octants — each face has exactly ONE outflowing group), where "reflect after the octant" == "reflect after the LAST octant outflowing the face". G-S gave the correct fixed point (rel-Linf 4.86e-8 vs Jacobi). The bug is INVISIBLE on axis-aligned quads — it requires a cubature where a face is shared by multiple octant groups.

2. **The discriminating signal was a closed-form anchor on a DIAGONAL cubature.** The all-reflective flat-flux anchor (`5.882`) is not invariant under the corrupted reflect, so it caught it (`3.4 ≠ 5.882`); a non-flat het case showed the opposite symptom (collapsed TO flat).

3. **The default flip exposed it.** The bug only reaches production when `inner_schedule="gauss_seidel"` is the 2-D default; the two lebedev tests had been green under the (pre-3c) Jacobi path. Wiring the selector + flipping the default surfaced both immediately — a worked instance of "a convention/default flip is the moment latent path-divergence declares itself".

**The fix.** Assign each reflective face to the LAST in-plane octant group (in sweep order) that outflows it (`last_group_for_face`), so the inter-group reflect fires only when the face's outflow is COMPLETE (every octant streaming out through it has been swept this pass). Octants reading the face that are swept BEFORE its reflect keep the lagged seed (the cyclic back-edge → partial one-pass G-S — always valid). Axis-aligned quads are unchanged (last == only). Restores the correct fixed point on lebedev (flat `5.882`; het ≡ Krylov).

**Which test catches it.** `tests/sn/sweep/core/test_sweep_schedule.py::test_gs_diagonal_quadrature_shared_face_assigned_to_last_group_only` (NEW, `@pytest.mark.catches("ERR-056")`) — a lebedev cubature schedule must assign each reflective face to EXACTLY ONE group, the LAST that outflows it (a direct structural pin). End-to-end: `test_2d_octant_sweep_closed_form_anchor` (lebedev all-reflective ≡ flat closed form) + `test_2d_heterogeneous_si_krylov_equivalence` (lebedev het SI ≡ Krylov) both gate the converged fixed point.

**Lesson.**

1. **A face shared by multiple work-units must be reduced only after the LAST contributing unit completes.** The boundary reflect is a REDUCTION over all octants outflowing a face; firing it after a partial contribution reads incomplete data. Whenever a schedule fans work over units that share an output slot, the slot's downstream consumer must be gated on the LAST writer, not the first. Generalises to any wavefront / fan-in scheduling (KBA, multigroup Gauss-Seidel over shared down-scatter targets, CP surface-current pairing).

2. **Ephemeral buffers make "stale" mean "garbage", not "last iterate".** A lagged (Jacobi) coupling reads the PREVIOUS iterate's real value, which is self-consistent at convergence. But a coupling that reads a REBUILT-each-solve buffer before the producer has run reads the SEED (an unrelated quantity) — and that error does NOT vanish at convergence. Distinguish "lagged but real" from "premature on ephemeral storage"; only the former is a valid splitting.

3. **A splitting/scheduling change MUST be gated on a structure that breaks the schedule's symmetry — here a DIAGONAL cubature.** Axis-aligned quadratures cannot see shared-face bugs (every face has one outflowing octant). The verification regime for any octant-schedule feature must include a cubature where octants outflow multiple faces (`lebedev` / `level_symmetric`), exactly as multi-group + heterogeneous is mandatory for scattering bugs. See `[[vv-principles]]` Mode 9 (synthetic-acceleration / splitting changes verified to the SAME fixed point on a stressing config).

**Test reference:** `tests/sn/sweep/core/test_sweep_schedule.py::test_gs_diagonal_quadrature_shared_face_assigned_to_last_group_only` (`@pytest.mark.catches("ERR-056")`); end-to-end gates `tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py::test_2d_octant_sweep_closed_form_anchor` + `tests/sn/solve/test_fixed_source_2d_equivalence.py::test_2d_heterogeneous_si_krylov_equivalence`. Fixed in `orpheus/sn/sweep_schedule.py` (commit `a39905a`).

## ERR-057 — `ordinate_scan` conditioning guard tested the `cumprod_a[-1] == 0` PROXY (exact reset) instead of the true non-finite condition, so a denormal-but-nonzero cumprod leaked a NaN through the Blelloch closed form

**Date:** 2026-06-10
**Module:** `sn` (`orpheus/sn/spatial/scan.py::ordinate_scan`; consumed by the 1-D sweep `_sweep_1d_unified` and — imminently — the S5 scan-march, one call per transverse line)
**Failure mode:** guard-predicate-incompleteness — a conditioning dispatch tested a PROXY for ONE cause of a failure instead of the failure CONDITION itself. Mechanically distinct from modes #1–#6: the code is not wrong on its own terms, the GUARD under-covers its own dispatch predicate.

**Symptom.** For a long, optically-thick chain with no exact reset (`a = 0.1`, `b = 1.0`, `nx = 312`, `|a| < 1` throughout), `ordinate_scan` returns NaN, while the mathematically-equivalent explicit serial fold stays finite (converging to the fixed point `b/(1−a) = 1.111…`). Demonstrated by the test-architect S5 probe: `cumprod_a[-1] = 1e-312` (a denormal, `!= 0.0`), `b / cumprod_a[-1] = inf`, `cumprod_a · inf = NaN`; the NaN then propagates wherever the scan output is consumed.

**Root cause.** `ordinate_scan` dispatches between the Blelloch closed form `cumprod_a · (psi_0 + cumsum(b / cumprod_a))` (fast, division-based) and the division-free pair-monoid backend (`_pair_monoid_scan`, exact + finite by construction). The guard was `if np.any(cumprod_a[-1] == 0.0)`, which tests for an EXACT reset (`a[i] = 0`, the cylindrical pole cell, ERR-054). But the closed form is ALSO invalid in a second regime the guard never covered: a long contractive chain decays `cumprod_a` into the IEEE-754 denormal band `(5e-324, 2.23e-308)` WITHOUT hitting exact zero, so `b / cumprod_a` overflows to `inf` and `cumprod_a · inf = NaN`. The `== 0.0` predicate is a PROXY for "the division will fail"; it under-covers, missing the denormal band where the division overflows but the divisor is not exactly zero.

**How it hid.**

1. **The exact-reset cause (ERR-054) was the only regime ever probed.** The reset suite drives `a[i] = 0` exactly (the pole-cell resonance), which the `== 0.0` guard DOES catch → pair-monoid → finite. The denormal band (`cumprod` small-but-nonzero) was never constructed in a test, so the gap was invisible.
2. **The fast-path bit-identity test uses short chains.** `test_fast_path_bit_identical` runs `nx ≤ 40` — far above the `cumprod[-1]` denormal threshold (`~0.1^308`), so the closed form stays well-conditioned and the gap never triggers.
3. **Production SN configs rarely reach `nx ~ 312` of near-constant attenuation.** Well-resolved DD keeps `|a| < 1` but typically not for ~300 cells of uniform thick attenuation on a single ordinate line — so the bug is reachable-in-principle but not yet observed in a production solve. The imminent S5 scan-march MULTIPLIES the exposure (one `ordinate_scan` call per transverse line × octant × group × iteration), which is why it was fixed BEFORE the carve built on top of it (Cardinal Rule 1, "fix bugs immediately").

**The fix.** Dispatch on the TRUE failure condition: compute the closed form (under an `np.errstate` that suppresses the expected intermediate inf/NaN), and fall to the pair-monoid iff `not np.all(np.isfinite(closed_form))`. This catches the exact reset AND the denormal underflow AND any cumsum overflow with one honest predicate, is bit-identical on every finite (i.e. all currently-passing) input, and costs one extra O(N) reduction over an array the closed form already materialised. `orpheus/sn/spatial/scan.py`.

**Which test catches it.** `tests/sn/spatial/test_ordinate_scan_reset.py::TestOrdinateScanDenormalUnderflow::test_denormal_cumprod_underflow_stays_finite` (`@pytest.mark.catches("ERR-057")`) — drives the chain into the denormal band (with a `pytest.fail` precondition asserting `cumprod[-1] ∈ (0, tiny)` so the regime cannot pass vacuously, `-O`-safe per Mode 8) and pins `ordinate_scan` finite + equal to the explicit serial loop. The bit-identity half stays pinned by `::TestOrdinateScanReset::test_fast_path_bit_identical`.

**Lesson.** A conditioning guard must test the failure CONDITION, not a PROXY for one of its causes. `cumprod == 0` is a proxy for "the division `b/cumprod` will produce a non-finite result"; the proxy under-covered the denormal band where the divisor is nonzero but the quotient still overflows. When a fast path is gated by a fallback condition, gate on the fast path's actual output validity (`np.all(np.isfinite(result))`), not on a hand-enumerated subset of the inputs that would invalidate it — the enumeration is an open set and will miss a case. Generalises to every backend-dispatch-by-conditioning (MOC optical-depth `exp` underflow, CP escape-probability small-τ series cutoff, any cumprod/cumsum closed form vs its stable rearrangement).

**Test reference:** `tests/sn/spatial/test_ordinate_scan_reset.py::TestOrdinateScanDenormalUnderflow::test_denormal_cumprod_underflow_stays_finite` (`@pytest.mark.catches("ERR-057")`); bit-identity guard `::TestOrdinateScanReset::test_fast_path_bit_identical`. Fixed in `orpheus/sn/spatial/scan.py` (issue #222, S5.0).


## ERR-058 — Curvilinear within-group closure SEEDS were self-referential / proxy-sourced: exact on flat ψ (every gate green), O(1)-wrong on every non-flat field — solution error floored mesh-independently

**Date:** 2026-06-12
**Module:** `sn` (`orpheus/sn/operator.py` `_compute_LpC` / `_compute_decomposition` / `_compute_LpC_transpose`; `orpheus/sn/loss_representation.py` curvilinear per-ordinate sweep; `orpheus/sn/spatial/psi_half_angle_seed.py`)
**Failure mode:** closure-seed inconsistency — two independent closure SEEDS in the curvilinear within-group operator were constructed so that they are EXACT on spatially/angularly flat fields and O(1)-wrong on every other field. Both pre-date the #195 investigation; both hid behind the flat-field exactness (vv-principles Mode 7 at the level of operator internals).

**Symptom.** The curvilinear isotropic L1 MMS L2 error PLATEAUS mesh-independently (sphere ~0.0413, cylinder ~0.0494, nx 20→640, orders → 0) on `main`, with SI ≡ Krylov bit-identical (post-unification: one discrete system, no independent inner to disagree). The pre-asymptotic-transient framing of issue #195 was empirically refuted: no refinement helps.

**Two manifestations (one class).**

**(a) Spatial pole-face seed (the +μ chain's r=0 inflow).** The matvec and sweep read the innermost CELL-CENTRE flux ψ(Δr/2) as the pole-FACE (r=0) value — a half-cell offset, O(h)-wrong on non-flat radial profiles, exact on flat ψ. The DD face chain propagates the seed error as an undamped odd–even alternation and the area-weighted streaming amplifies it by ~A/V. **Fix:** the Carlson coupled-pole seed — at r=0 the outward characteristic CONTINUES the inward one, ψ(0, +μ) = ψ(0, −μ), so the −1 sweep runs FIRST and its pole-face outflow at the mirror ordinate (`quad.reflection_index("x")`) seeds the +1 sweep. Data, not self-reference; lower-triangular; flat-flux exact; the #192-deferred "inward-determines-outward" pole condition. Adjoint: the +1 seed cotangent routes into the −1 reversal's initial outflow cotangent at the mirrored ordinates.

**(b) Angular half-angle thread seed (the M-M recurrence's φ̄₁/₂).** `CarlsonInwardSweep` solved the Hébert §3.9.4 starting-direction ODE with the PROXY source Q̄ = Σ_t·φ₀/Σw — equal to the true within-group source ONLY at flat-flux equilibrium (Σ_t φ₀ = Q̄). On any non-equilibrium field the seed solves the wrong ODE: measured φ̄ = 0.5777 where the isotropic input value is 0.5000; per-ordinate redistribution residual up to ±55 at the pole, ±13 in the bulk (vs continuous ±0.31). **This was the dominant defect.** **Fix:** `AngularEdgeExtrapolation` — the half-angle seed is, definitionally, the INPUT FIELD's value at the level's angular edge μ_start (sphere −1; cylinder −√(1−ξ²), single-sourced from the α/τ construction site in `reduced_operator.py` and threaded via `CarlsonSweepContext.mu_start`, a REQUIRED field). Two-point linear extrapolation in μ through the level's two most-inward distinct-μ ordinates: exact on angle-flat fields, exact on linear-in-μ fields (and the M-M recurrence with unclamped τ then threads the whole half-angle grid exactly), O(Δμ²) generally, linear in the input, no σ_t/source/BC dependence. The seed strategy now OWNS its adjoint (`seed_adjoint`) so `angular_adjoint` delegates instead of hardcoding a twin.

**How it hid (both manifestations).**

1. **Flat-field exactness.** Cell-centre = face value on flat ψ; Σ_tφ₀ = true source at flat-flux equilibrium. Every per-ordinate flat-flux identity gate, every streaming-equilibrium L0 gate, and every reflective-equilibrium configuration sits EXACTLY in the blind regime — Mode 7 operating on operator internals rather than on an MMS ansatz.
2. **α-dome telescoping made scalar checks blind to (b).** Σ_n w_n of the M-M redistribution telescopes to zero REGARDLESS of the half-angle values, so any weight-summed (scalar/balance) residual cannot see a wrong half-angle thread — anti-pattern #8 instantiated inside a diagnostic: the #195 probe-3 scalar residual went O(h²) after fixing (a) alone while the per-ordinate residual was still O(10), which mis-supported a "near-singular operator" hypothesis until the per-ordinate check (σ_min ≈ 0.9 — never near-singular).
3. **Historical compensation.** The Phase-D-era apply path measured O(h²) under Krylov (issue #195's premise) — its closure internals compensated differently; the Depth-B/Wave-T matvec rebuild (window `c4fa030..90e7d4e`) changed the redistribution assembly and surfaced the latent seed inconsistency. SI (the sweep) had ALWAYS plateaued (#195's own data) — the wrong fixed point was the SWEEP's all along (#98's original 35%-at-r=0 finding is this same class).

**Post-fix evidence (2026-06-12).** Sphere isotropic MMS: errors [1.49e-2, 3.73e-3, 9.28e-4, 2.31e-4, 5.74e-5] at nx ∈ {20..320}, orders 2.00–2.01; cylinder [2.16e-3 … 3.37e-5], orders 2.00; SI ≡ Krylov bit-identical; magnitude band `1e-8 < err < 1e-3` satisfied (sphere nx ≥ 80, cylinder nx ≥ 40). Per-ordinate volume-weighted operator residual on ψ_ref decays (sphere order ~1.5 — the pole-adjacent bounded band under the r²dr volume weight; cylinder order 2.0, pointwise O(h²) everywhere). Anisotropic cases dropped ~50× and are now limited by the fixed-quadrature angular floor of the half-angle thread interpolation — a test-design retune tracked at #229, NOT a residual instance of this class.

**Which tests catch it.** `tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py::test_operator_admits_isotropic_mms_per_ordinate` (`@pytest.mark.catches("ERR-058")` — the fast per-ordinate volume-weighted admission gate, the structurally decisive check); `tests/sn/verification/mms/test_mms_curvilinear.py::test_sn_{spherical,cylindrical}_mms_converges_second_order` (`catches("ERR-058")` — the end-to-end L1 ladders whose xfail markers came off with this fix); the flat-flux and streaming-equilibrium gates pin the flat-field exactness both fixes preserve; `tests/sn/operators/test_g_adjoint_reciprocity.py` pins the strategy-owned seed adjoints.

**Lesson.** A discrete CLOSURE SEED is part of the operator and must be verified per-ordinate on a NON-FLAT field: (i) any seed that is a self-referential read of the unknown family it seeds (cell-centre-as-face) or a proxy-sourced auxiliary solve (Σ_tφ₀ for the true source) should be treated as wrong-until-proven on non-equilibrium data; (ii) scalar/balance residuals are structurally blind to redistribution-thread defects (the telescoping is BY CONSTRUCTION) — residual audits of angular closures MUST be per-ordinate; (iii) "every gate is green" is evidence about the gates' regime, not about the operator: enumerate which fields each closure seed is exact on, and place at least one gate outside that set.

**Test reference:** `tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py` (`catches("ERR-058")`), `tests/sn/verification/mms/test_mms_curvilinear.py` (both ladders, `catches("ERR-058")`). Fixed in `orpheus/sn/operator.py`, `orpheus/sn/loss_representation.py`, `orpheus/sn/spatial/psi_half_angle_seed.py` (+ `mu_start` threading through `orpheus/geometry/reduced_operator.py` → `StreamingTerms` → `GeometryCoefficients` / closure binding). Issues #195/#168/#99/#98; follow-up #229.

---

## ERR-059 — Curvilinear central-cell (r→0) scalar flux is first-order O(h): inherent to plain diamond difference, invisible to the volume-weighted L2 / k_eff gates (documented limitation, WONTFIX-for-DD)

**Date:** 2026-06-13 (W2 of the curvilinear-aniso program, Issue #233)
**Module:** `sn` (the curvilinear diamond-difference cell update — `orpheus/sn/spatial/diamond.py` / `cell_balance.py` central-cell branch; the `dd-curvilinear-scalar` cell-update closure).
**Failure mode:** #5 Index/closure error at a coordinate singularity — the diamond closure `ψ̄ = ½(ψ_in + ψ_out)` is inconsistent at the zero-area inner face A(0) = 0, giving first-order spatial accuracy at the single central cell. NOT a wrong-answer-accepted bug (ERR-026 was that, now CLOSED) — this is an INHERENT property of the standard scheme.

**The defect.** At the r → 0 central cell the inner radial face has zero area (A(0) = 0). The diamond spatial closure then gives `ψ_out = 2ψ̄`, which **over-predicts the pole outer face by exactly +50 %** (mesh-independent rel. error 0.5000) relative to the true face value `A(h)` vs the volume-average `2⟨A⟩_vol = 2·¾A(h) = 1.5·A(h)`. Deeper: the conservative *balance itself* is inconsistent at the pole — fed the EXACT cell average and EXACT inflow it solves for an outer face −46 % wrong, and the residual-per-volume plateaus mesh-independently, because `A_in = 0` degenerates the streaming surface integral while `V ~ h³`. Result: the central-cell scalar flux converges at O(h^0.95) (pointwise / L∞), not O(h²), while the interior converges clean O(h²).

**Decomposition (W2 numerics-investigator, diag_16..31).** The apparent error is ~75 % a comparison artifact + ~25 % genuine inherent residual:

1. **~75 % MMS comparison artifact.** The production spherical MMS evaluates the source at the cell MIDPOINT and compares `phi_solver` against `phi_exact(midpoint)`. But the spherical DD discrete unknown IS the cell-VOLUME average `(4π/V)∫r²φ dr` (Hébert 2009 Eq. 3.430 — the unknown is *defined* as the shell average, not a point value). Under r²dr weighting the volume-average and the midpoint differ by O(h) at the pole (volume-centroid at ¾h, not ½h). Using the shell-averaged source AND the shell-volume-average reference drops the pole error ~4× (0.0212 → 0.00497).
2. **~25 % genuine but LITERATURE-ACCEPTED INHERENT first order.** Even the fully-consistent FV-MMS (shell-avg source + shell-avg ref) leaves the pole at clean O(h^1.00) — the diamond-at-zero-inner-area inconsistency above. Hébert §3.9.4 and Stacey §9.9 BOTH use this plain-diamond + Carlson-starting-direction + symmetry scheme at the central cell with NO special O(h²) closure and NEITHER flags reduced order there (literature-researcher: Tier-2 search found no paper on a special O(h²) central-cell closure).
3. **NOT cleanly fixable by a local closure.** The volume-weighted linear reconstruction `ψ̄ = β·ψ_out + (1−β)·ψ_in` with β = ¾ at the pole (the value O(h³)-consistent vs ⟨A⟩_vol at exact faces) was validated end-to-end (faithful production-sweep monkeypatch + β=½-identity guard to 3e-16): β = ¾ does NOT restore O(h²) — pole stays O(h), magnitude slightly worse, full-mesh β degrades the interior. A genuine fix needs a non-local higher-order central-cell scheme: linear-discontinuous (#6), cell-updates (#158), or nodal (Wu-Xie-Fischer 1999 NSE 133).

**Why it hid (and why it needed a new probe to surface).** It is INVISIBLE to the production volume-weighted L2 norm AND to k_eff:

- The production `test_sn_spherical_mms_converges_second_order` uses the volume-weighted L2 `√Σ V·diff²`. The pole O(h) at ONE cell of V ~ h³ contributes √V ~ h^1.5 → ~h^2.5 to L2 — subdominant to the interior O(h²). Both midpoint AND volume-average L2 references converge clean O(h^2.00); only L∞ (pole) is O(h).
- k_eff is shape-blind: reflective sphere = k_inf = 1.875 EXACT mesh-independent; vacuum sphere converges monotone to ~1.78590 at O(h^1.48), increments 2e-5 at nx=160 (far below engineering tolerance).
- This is the SECOND time a curvilinear defect hid behind a norm choice: ERR-058 hid behind flat-field exactness; ERR-059 hides behind the L2 volume weight + the midpoint comparison. The decisive surfacing probe was a per-cell pointwise/L∞ rate decomposition (`E_test = E_artifact(midpoint−volavg) + E_true(solver−volavg)`).

**Cylinder shares the same defect, MASKED.** Cyl pole vs MIDPOINT is O(h²) (1.94/1.97/1.98) but vs the VOLUME-AVERAGE is O(h) (0.99/0.99/1.00) — the SAME diamond inconsistency, accidentally hidden by the midpoint comparison (the cyl r dr linear weight puts the volume-centroid where diamond's ½A(h) ≈ midpoint A(h/2)). So the cylinder pole is NOT "clean O(h²)" — it is the same O(h) volume-average defect.

**Which test catches it.** `tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py` (W2, commit `255eba4`) — a CHARACTERIZATION gate (`@pytest.mark.l1`, `@pytest.mark.catches("ERR-059")`). It pins the limitation WITHOUT calcifying it: the pole L∞ order is **lower-bounded only** (> 0.8, "at least first order, does not regress"), the pole is the L∞-dominant cell (fraction > 0.99), and the interior is clean O(h²) (> 1.8) — NO upper bound on the pole order, so a future LD/nodal fix that lifts the pole to O(h²) keeps the gate green (vv anti-patterns #5/#17). The GUARANTEE half (carrying `verifies("dd-curvilinear-scalar", ...)`) asserts the global volume-weighted L2 is O(h²) under BOTH a midpoint AND a Hébert-3.430 shell-volume-average reference (built from `scipy.integrate.quad` — a trusted-library integrator structurally independent of the solver); agreement on the order across two structurally-different references proves the L2 order is REAL, not a midpoint artifact.

**Status: DOCUMENTED INHERENT LIMITATION 2026-06-13 (WONTFIX-for-DD).** Issue #233 stays **OPEN** to track the future higher-order central-cell scheme (LD #6 / cell-updates #158 / nodal). NOT a code bug in the diamond-difference scheme; it is the accepted, literature-unflagged behaviour of plain diamond at a coordinate singularity. Distinct from #168 (outer face, CLOSED), ERR-058 / #195 (the closure seed, CLOSED), and the angular floors (#229). Full theory treatment: `docs/theory/methods/sn/index.rst` §`sn-pole-cell-spatial-closure`.

**Lesson.** A convergence-rate gate on a volume-weighted L2 norm is structurally blind to an error concentrated at a measure-zero-in-the-limit cell (the pole, V ~ h³). When verifying a curvilinear / coordinate-singularity scheme, ADD a pointwise / L∞ per-cell probe alongside the integral-norm gate — and choose the MMS reference quantity to MATCH the scheme's discrete unknown (the cell-VOLUME average per Hébert 3.430, NOT a midpoint point-value), or the comparison itself injects a spurious O(h) at the singular cell. The shell-volume-average reference is the structurally-independent ground that separates the comparison artifact from the genuine residual.

---

## τ-clamp mis-citation finding (W1, 2026-06-13) — a correctness-of-design finding, not a bug-class instance

**Date:** 2026-06-13 (W1 of the curvilinear-aniso program; companion to ERR-059).
**Module:** `sn` / `geometry` (`orpheus/geometry/reduced_operator.py` `spherical_streaming`).
**Class:** correctness-of-DESIGN finding (a mis-cited, over-conservative patch that was 100 % spurious on physical fields) — NOT a new wrong-answer instance, so it carries no `catches("ERR-NNN")` of its own. It is the W1 leg of the ERR-026-family clean-up.

**Finding.** The spherical Morel--Montry weighted-diamond weight `τ_n = (μ_n − μ_{n-1/2})/(μ_{n+1/2} − μ_{n-1/2})` had been wrapped in a `[½, 1]` clamp, cited to Lewis & Miller §4.5. Triple-confirmed mis-citation:

1. **Literature.** Bailey-Morel-Chang (2010) NSE 165 **Eq. 43** give exactly this UNCLAMPED `τ_raw` as the unique weight exact for a flux linear in μ, admissible range τ ∈ [0, 1]; Hébert §3.9.4 uses pure diamond (τ = ½), no clamp; Lewis & Miller §4.5 does NOT prescribe the `[½, 1]` clamp — the citation was wrong.
2. **Positivity never needed.** On every realistic converged solve (smooth MMS, homogeneous k_eff = 1, thick absorber) there are ZERO negative half-angle fluxes, clamped or unclamped; every clamp activation is spurious (160/320/80/240 activations across stress configs, 0 protective). On GL quadrature `τ_raw ∈ [0.39, 0.61]` (never 0), always interior to [0,1].
3. **Stability without it.** Unclamped sphere SI converges with strictly positive, finite scalar flux on every stress config (thick absorber, near-vacuum, c = 0.999, S64).

**Fix (W1, commit `b2d8a6d`).** Drop the clamp for the SPHERE only (use raw `τ_raw`); CYLINDER keeps it (product/LS quadratures put the most-inward azimuthal ordinate exactly on η = −sinθ → `τ_raw = 0` exactly, a structural ÷0 block needing a 2-D (η,φ) closure, #229). The removal is STATIC (config-time), so the operator stays linear and the SI ≡ Krylov twin identity holds — a *dynamic* (ψ-dependent) negative-flux fixup would make the operator nonlinear and break the twin.

**Mixed accuracy signature (the gotcha worth recording).** Unclamping cleans the coarse rate (sphere S16 coarse orders 1.978 → 1.995) but RAISES the S16 fine floor (7.3e-4 → 1.2e-3): the lower clamped floor was a *fortuitous cancellation* (the clamp's constant bias partly offset the #229 angular-thread interpolation floor), not a genuine accuracy gain. Iso solves are unchanged in real arithmetic but NOT bit-identical at IEEE-754 (the τ-bearing reduction returns ψ exactly only ~81 % of the time, ≤1 ULP else; converged drift |Δk| = 2.3e-13 on homogeneous-reflective sphere — an FP-tail anchored to k_inf = 1.875).

**Which gate.** `tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py` (closure-unit τ-independence on flat fields; converged iso anchored to k_inf; unclamped positivity) + the W1 `@slow` aniso discriminators in `tests/sn/verification/mms/test_curvilinear_aniso_convergence.py`.

**Lesson.** A positivity clamp on a discretisation weight is a tax, not a safeguard, UNTIL you have measured that the negativity it guards against actually arises on physical (converged, smooth) fields — and traced the clamp to a primary source. A clamp justified only by "it keeps the closure positive" against fields the physics never produces should be treated as wrong-until-the-negativity-is-demonstrated. (Companion to ERR-058's "every gate green is evidence about the gates' regime, not the operator" — here: every clamp activation is evidence about pathological roughness, not physical solves.)

## ERR-060 — d-generic UBLD symbolic inflow assembler dropped the |μ_axis| streaming factor: passed the d=1 reduction oracle, FAILED exact-on-bilinear (Oracle ii)

**Date:** 2026-06-16 (D5b-S1 Branch 1, the d-generic UBLD algebra-of-record; #240/#38/#37).
**Module:** `sn` derivations (`orpheus/derivations/discrete/sn/ld_ubld.py` — `assemble_inflow_axis`).
**Class:** caught-at-derivation (Mode-3 missing factor) — the bug never reached production; the symbolic oracle caught it.

**Failure mode:** #3 missing factor. The first draft of `assemble_inflow_axis` dropped the `|μ_axis|` streaming factor on the inflow RHS of the d-generic UBLD cell system. The d=1 reduction oracles build the d=1 inflow RHS INLINE (not through `assemble_inflow_axis`), so they were structurally blind to the dropped factor.

**Impact.** None in production (caught in the Branch-1 symbolic reference before any numpy code consumed the assembler). Had it shipped, the d≥2 UBLD cell inflow would have been under-weighted by `1/|μ_axis|` per axis.

**How it hid.** All three d=1 oracles passed (d=1 RHS built inline). The d=2 exact-on-bilinear gate (Oracle ii) is the FIRST consumer of the multi-axis inflow assembler — it manufactured a known bilinear `a+bx+cy+dxy`, derived `Q = Ω·∇ψ + Σ_t ψ` via SymPy, solved the 4 cell moments, compared to the exact bilinear; the dropped `|μ_axis|` made the solved moments differ (the xy coupling is active, `d≠0`).

**Fix.** `return mu_axis * sp.Matrix(...)`. Mutation-verified `-O`-safe: re-dropping the factor makes `test_d2_exact_on_bilinear` FAIL (returncode 1, via `pytest.fail`) while the d=1 oracles stay GREEN.

**Which test catches it.** `tests/sn/spatial/test_ld_ubld_symbolic.py::test_d2_exact_on_bilinear` (Branch 1) + `tests/sn/spatial/test_ld_ubld_primitive.py::test_d2_exact_on_bilinear` (Branch 2 numpy), both `@pytest.mark.foundation @pytest.mark.catches("ERR-060")`, Mode-8-safe (`pytest.fail`). NOTE: `tests/sn/spatial/test_linear_discontinuous.py::test_d2_assembled_matrices_match_symbolic` carried `catches("ERR-060")` but is BLIND to it (it checks `assemble_ubld`'s A/M/G/F_out, which carry no inflow factor, and PASSES under the |μ_axis| drop). Only the exact-on-bilinear gates are genuine catchers (the marker on the A==A pin is a coverage-claim error to drop — qa L-031).

**Lesson.** A reduction oracle that builds its own reduced-case RHS inline is blind to the general assembler's higher-d terms — ship the exact-on-bilinear (d≥2) gate WITH the d=1 reduction oracle, never the d=1 reduction alone. The d=1-blind / d≥2-caught split is the H2 signature lifted to dimension. → numerical-bug-signatures Signature 4 family at the symbolic-derivation layer.

## ERR-061 — LD slope moment stored in the per-ordinate SWEEP frame, not the GLOBAL-x frame: the angular reduction φ̂ = Σ_n w_n ψ̂_n cancels forward against backward ordinates → the diffusion-limit slope source Σ_s·φ̂ is ~6× under-driven → LD loses the thick-cell diffusion limit

**Date:** 2026-06-17 (D5b-S3, the Increment-C moment-iterate fold; #240/#38/#37).
**Module:** `sn` (`orpheus/sn/sweep_graph.py` `_CellSolve` / `_CellResidual`; `orpheus/sn/loss_representation.py` the d=1 matvec `_apply_walk`).
**Class:** wrong-answer — the operator was internally self-consistent (matvec≡sweep round-trip 1e-16, SI≡Krylov on the SAME operator) but converged to a DIFFUSION-INCONSISTENT fixed point.

**Failure mode:** #1 sign flip + #6 convention drift. The per-cell LD kernel produces/consumes the `2^d` moment vector in the per-ordinate SWEEP frame (each axis oriented so the downstream face is at local `+1`). For an ordinate sweeping in the NEGATIVE global direction on axis `a` (`octant_signs[a] == -1`) the sweep coordinate is the REVERSE of the global coordinate, so the slope (`P₁`) moment on that axis is sign-FLIPPED relative to the global-x slope. The iterate `φ̂` + its isotropic scattering source `Σ_s·φ̂` live in the GLOBAL frame: `φ̂ = Σ_n w_n ψ̂_n` sums slopes across ordinates of BOTH sweep directions. The producer (`_CellSolve.cell` emit) stored the raw sweep-frame slope and the consumer (`integrate_angular` / `S.apply`) summed it as if global-frame — so the backward ordinates' opposite-signed slopes partially CANCELLED the forward ones.

**Impact.** On a thick scattering slab (σ_t=40, c=0.99, σ_t·h=10/cell, vacuum) at nx=4, LD gave 1.4717 vs DD 2.4069 (rel 38.9%) and vs the analytic diffusion solution 2.31 (rel 36%). The converged scalar slope `φ̂` was ~6× too small to satisfy the LM-1989 discrete diffusion continuity (Eq 4.9b, `φ̄_j + φ̂_j = φ̄_{j+1} − φ̂_{j+1}`); the per-ordinate `ψ̂_n` at a cell with a positive global-x gradient had `+0.048` (forward) but `−0.028` (backward) — opposite signs, the smoking gun. The error did NOT grow with refinement (it shrank as cells thinned and the slope source mattered less) — the classic flat-source-LD signature persisting THROUGH the slope-source thread.

**How it hid (the instructive part).** Every component was individually CORRECT: the per-cell LD 2×2 matched LM-1989 Eq (4.3) exactly; the dense UBLD solve matched the analytic 2×2; the scattering produced `Σ_s·φ̂/W` at full strength on the slope row (ratio 1.0); the matvec round-trip vanished to 1e-16; SI≡Krylov on the SAME operator (D5b-S3 GATE 3, the genuine Mode-9). The bug was the FRAME-CONSISTENCY between two correct components — a wrong-fixed-point that the matvec-self-consistency gate is structurally BLIND to (it proves SI and Krylov solve the SAME operator; it CANNOT prove that operator is diffusion-consistent — vv-principles §5, "O(h²) to the wrong limit is still O(h²)"). The decisive evidence was a STRUCTURALLY-INDEPENDENT from-scratch LD-SN solver (a direct LM-1989 2×2 + source iteration, NO ORPHEUS kernel) that reproduced ORPHEUS's WRONG value bit-for-bit when it summed sweep-frame slopes, and RECOVERED the diffusion limit (nx=4 rel 2.3%) when it stored global-frame slopes — pinning the root cause to the slope-moment frame independent of ORPHEUS's code.

**Fix.** A single-sourced `2^d` moment-frame involution `octant_moment_frame_signs(octant_signs, per_axis)` = `∏_a (octant_sign_a)^{o_a}` (average moment sign-invariant; per-axis slope flips once if that axis sweeps backward; the d=2 cross moment `x̂y` flips when an ODD number of its axes reverse). Applied via the `_reframe` helper at both cell ops: the source/probe is mapped global→sweep on INPUT and the emitted moment/residual sweep→global on OUTPUT (the map is its own inverse). The OUTGOING FACE (`psi_out`) stays sweep-frame — it propagates along the wavefront and never crosses into the global-frame iterate. DD/Step (`per_axis == 1` → `None`) are byte-identical (the negative control: GATE 4 = 513 pass / 1 skip / 4 xfail, zero drift). The flat scalar source (matvec zero / flat external — only the sign-invariant average moment) is frame-invariant and skipped by the `arr.shape[-1] != frame_signs.shape[0]` guard, so it is never broadcast into a spurious moment axis. Post-fix: nx=4 LD vs DD rel 38.9% → 4.1%; nx=16 7.9% → 0.2%; nx=64 0.9% → 0.0%. The 2-D analog converges 8.4% → 1.7% → 0.4% across n=4/8/16.

**Which test catches it.** `tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit` (1G) + `::test_ld_thick_diffusive_limit_2g` (2G-het, Mode-6 group-coupled slope source) — both `@pytest.mark.l1 @pytest.mark.catches("ERR-061")`, both Mode-8-safe (`np.testing.assert_array_less`, fires under `-O`). The slope-frame fingerprint is pinned by `derivations/diagnostics/diag_240_d5b_s3_probe_11_root_cause.py` (forward and backward ordinate slopes must share sign in the global frame) and the structurally-independent confirmation by `diag_240_d5b_s3_probe_08_independent_ld.py` (from-scratch LD recovers diffusion only with the global-frame correction).

**Lesson.** A per-ordinate spatial-moment quantity (an LD slope, an Pℓ anisotropic moment) that is produced in a direction-dependent SWEEP frame MUST be lifted to the global frame BEFORE the angular reduction that sums it across ordinates — the producer and the consumer must agree on the frame, or forward and backward ordinates cancel a quantity that should reinforce. The matvec-self-consistency gate (SI≡Krylov, round-trip≈0) is necessary but NEVER sufficient for a moment-iterate fold: it proves the operator is internally consistent, not that its fixed point is the physically correct one — gate the converged VALUE against a structurally-independent reference (here: the continuous diffusion solution + an independent from-scratch LD kernel), never the round-trip. → numerical-bug-signatures: a NEW frame-convention class (sweep-frame vs global-frame for direction-dependent moments) adjacent to Signature 3 (scattering transpose) and Signature 4 (quadrature normalization) — the common thread is a per-ordinate convention that is invisible until a quantity is summed across ordinates of opposite sweep direction (the angular reduction is the discriminator, exactly as H2/H3 predict: flat flux nulls the slope, and conservation/round-trip are telescoping-degenerate to the frame error).

## ERR-062 — LD matvec pure-z arm `(L+C)ψ̄ = σ·ψ̄` lacked the moment-axis broadcast guard its sweep twin `ψ̄ = Q/σ_t` already had → `inner_solver="krylov"` on a 2-D LD mesh with a pure-z-bearing quadrature broadcast-crashed

**Date:** 2026-06-17 (D5b-S3, closing the two review findings before commit; #240/#38/#37).
**Module:** `sn` (`orpheus/sn/loss_representation.py` — the `loss_action` matvec `pure_z` arm vs the `octant_jacobi` sweep `pure_z` arm).
**Class:** crash (Mode-2 / Mode-3) — `ValueError: operands could not be broadcast together with shapes (ng,nx,ny) (1,ng,nx,ny,2^d)`; not a wrong-answer, a hard failure on a live config.

**Failure mode:** #2 variable/path swap (the L21 twin-path asymmetry) + #3 missing factor (the moment-broadcast reshape, present in one twin, absent in the other). The pure-z degenerate ordinates (μ_x = μ_y = 0, the ±z poles) have no in-plane streaming, so the cell is collision-only: the SOLVE balance is `ψ̄ = Q/σ_t` and its matvec twin is `(L+C)ψ̄ = σ_t·ψ̄`. At a multi-moment closure (LD, #240 D5b-S3) the source `Q` / probe `ψ̄` carry a trailing `2^d` spatial-moment axis that `σ_t` `(ng, *spatial)` lacks, so `σ_t` must gain a length-1 trailing axis to broadcast. The SWEEP arm HAD the guard (`if q.ndim > sig.ndim + 1: sig = sig[..., None]`); the MATVEC arm wrote the bare `sigma * probe[oct_idx]` — so a moment-valued probe broadcast-FAILED.

**Impact.** `solve_sn_fixed_source(scheme=LinearDiscontinuous(), inner_solver="krylov")` on ANY 2-D Cartesian LD mesh whose quadrature carries pure-z ordinates raised `ValueError` at the FIRST Krylov matvec. The production 2-D Cartesian MMS uses `Quadrature.lebedev(order=17)` (N=110) which has exactly 2 pure-z ordinates; `level_symmetric` has NONE, which is why the bug hid through the whole D5b-S3 development (every committed 2-D LD test used `level_symmetric`). The SI path was unaffected (its sweep arm had the guard); only the Krylov/matvec path crashed.

**How it hid (the instructive part).** It is the canonical L21 twin-path asymmetry: the sweep `Q/σ_t` and the matvec `σ_t·ψ̄` are two applications of the SAME collision-only operator, but they were written as separate `pure_z` closures, so the moment-broadcast convention was added to one and forgotten on the other when the `2^d` axis landed (D5b-S3 OWED-2 / the unified matvec). The matvec needed a COMMITTED gate, not a round-trip — the round-trip and the FFW≡MFW gates ran on `level_symmetric` and never exercised the pure-z arm at all. (Recurring a THIRD time per the L14/L18/L21 note: "the matvec needs a committed gate.")

**Fix.** Single-source the moment-broadcast through `_moment_broadcast_sigma(sig, moment_valued)` (= `sig[..., None] if moment_valued.ndim > sig.ndim + 1 else sig`), called by BOTH the sweep arm (`Q / _moment_broadcast_sigma(sig_t, q)`) and the matvec arm (`_moment_broadcast_sigma(sigma, probe_oct) * probe_oct`). The two twins now CANNOT diverge on the moment-axis reshape (Pattern 2 / L21 — sweep ≡ matvec are two applications of the same collision-only operator). DD/Step (no moment axis) → `sig` unchanged, byte-identical (the negative control: GATE 4 = 513 pass / 1 skip / 4 xfail, zero drift).

**Which test catches it.** `tests/sn/verification/mms/test_mms_ld_2d.py::test_ld_2d_krylov_equals_si_pure_z_quadrature` (`@pytest.mark.foundation @pytest.mark.catches("ERR-062")`, Mode-8-safe — `np.testing.assert_allclose` / `pytest.fail`). Mode-9 degeneracy-break config: a pure-z-bearing **Lebedev order 5** quadrature (N=14, genuine `mu_y` + the 2 ±z poles), **heterogeneous** 2-material map (A/B), **2-group asymmetric** XS with **non-zero self-scatter**, **NON-SQUARE** 5×4, vacuum edges. Mutation-verified: re-introducing the bare `sigma * probe[oct_idx]` makes the gate FAIL with the exact `ValueError` (shapes `(2,5,4)` vs `(1,2,5,4,4)`); with the fix Krylov ≡ SI to ~1e-11 (the same `(L+C−S_full)` fixed point).

**Lesson.** When a per-cell quantity gains a new axis (here the `2^d` spatial-moment axis), every TWIN that touches it — the sweep solve AND its matvec apply — must learn the new broadcast convention at the SAME single-source site, or the twin that was edited diverges from the twin that was forgotten. The L21 rule "sweep and matvec are two applications of the same operator" is not just a correctness aesthetic: a collision-only `Q/σ_t` ↔ `σ_t·ψ̄` pair that does not share its σ-reshape helper is a crash waiting for the first quadrature that exercises the degenerate arm. And the matvec ALWAYS needs a committed gate on the bug's exact habitat (a pure-z-bearing quadrature here), never a round-trip on a degenerate config — the round-trip is structurally blind to an arm it never runs. → coding-elegance Pattern 2 (single source) + the recurring L14/L18/L21 "the matvec needs a committed gate."

## ERR-063 — "zeroing χ on non-fissile materials is inert" assumed the SN/`compute_macro_xs` fission contract (χ gated by the SAME region's νΣf) and was FALSE for `solve_peierls_mg`, whose fission operator weights source-region νΣf by the SINK-region χ → zeroing region-B/C/D χ silently changed 7 L1 Class-B heterogeneous k_eff pins

**Date:** 2026-06-21 (#257 S10a QA validation; branch `feature/field-typed-operator-algebra`, HEAD `c6e21c0`, NOT committed).
**Module:** `data` + `derivations` (the `xs_library.py` non-fissile-χ precursor vs the `peierls_nystrom/geometry.py::solve_peierls_mg` MG fission operator).
**Class:** silent wrong-answer (Mode-3 missing-factor analog) — a behavior-neutrality premise that held for ONE fission contract and failed for ANOTHER consumer of the same `chi` field.

**Failure mode:** an unverified behavior-neutrality assumption. The S10a precursor zeroed the (placeholder, non-zero) χ on the NON-fissile xs_library regions B/C/D so the new `Mixture/Isotope.__post_init__` `assert_null()` guard would not reject them. The stated justification: "the fission source `χ·(νΣf·φ)` is zero regardless of χ because νΣf≡0 on a non-fissile material." TRUE for the SN `FissionOperator` and `compute_macro_xs` (χ is applied to the SAME material's νΣf). FALSE for `solve_peierls_mg`, whose multi-region fission operator is `B[i·ng+ge, j·ng+gs] += K[i,j] · chi[i,ge] · nu_sigma_f[j,gs]` — the emission spectrum is the SINK node `i`'s χ, the production is the SOURCE node `j`'s νΣf. A fission neutron born in fissile region A (node `j`) and emitted INTO non-fissile region B (node `i`) is weighted by χ_B. Zeroing χ_B removes that emission path, so k_eff moves O(1).

**Impact.** Direct behavior-neutrality probe on `solve_peierls_mg` (2R cylinder, `white_hebert`, `n_bc_modes=1`): region-B χ `[1,0]→[0,0]` moves k_eff `1.0985→0.5563` (1G/2R) and `1.1008→0.3856` (2G/2R). 7 L1 tests in `tests/derivations/test_peierls_rank_n_class_b_mr_mg.py` (cylinder hebert overshoot 1G/2R + 2G/2R; sphere hebert overshoot; `recovers_kinf[2G_2R]` + RICH; `2g_2r_rank1_mark_floor_pinned[cylinder-1d]` + `[sphere-1d]`) FAIL under S10a (errors swing to −44%…−75%, opposite sign from the expected +10%/+45% overshoot pins) but PASS at clean HEAD (proven via a `git worktree add c6e21c0` run: 4 + 4 passed). Only the RICH variant is `@slow`; the other 6 are plain `@pytest.mark.l1` and run in any `not slow` L1 suite. NOT a guard ValueError — a silent numerical regression.

**How it hid (the instructive part).** The closeout's behavior-neutrality evidence was (a) the DD-regression snapshots "did not move" and (b) a full-suite run with "0 EmissionSpectrum reds." Both are blind to this bug: (a) the DD regression only builds SN solves (the SN fission contract IS χ-gated-by-same-region, genuinely inert — my direct byte-identical re-solve confirmed keff `1.2298233055738448` bit-for-bit with χ_B = 0 vs [1,0]); (b) "0 EmissionSpectrum reds" counts only `assert_null`/`assert_normalized` ValueErrors — a silent accuracy-band regression in a DIFFERENT solver raises no EmissionSpectrum error. The slow Class-B peierls suite (494 s) was never run. The xs_library region dicts are a SHARED source consumed by BOTH the χ-gated SN path AND the sink-region-χ peierls path; a behavior-neutrality claim proven for one consumer was assumed for all. (The test authors had already flagged this exact χ-dependence in commit `76b11e8`: "chi-dependence of Hébert error reframes 2G/2R near-exact claim, Issue #132".)

**Fix (recommended; not applied at QA time).** Do NOT zero the shared xs_library B/C/D χ. The χ value on a non-fissile region is a LEGITIMATE input to the sink-region peierls fission operator. Options: (i) give the peierls MG cases their own χ source decoupled from the guarded library, or fix `solve_peierls_mg` to use the SOURCE region's χ (the physically-correct emission convention — a separate correctness question, file an issue); (ii) make the S10a guard key on production (`SigP` / `nu*sig_f`) rather than `SigF`, AND only apply `assert_null` to fields that are genuinely consumed as same-region χ·νΣf; (iii) narrowest: leave the library χ intact and restrict the guard so it does not force `assert_null` on placeholder library regions. The S10a TYPE + guard + SN-path behavior-neutrality are sound; only the "zero the shared library χ" precursor is unsafe.

**Which test catches it.** Existing `tests/derivations/test_peierls_rank_n_class_b_mr_mg.py::test_class_b_cylinder_hebert_heterogeneous_overshoot_known` (+ the 6 sibling Class-B heterogeneous pins) — they go red under the precursor and green at clean HEAD. The decisive structurally-independent probe is the direct `solve_peierls_mg` k_eff comparison (region-B χ `0` vs `[1,0]`), which shows an O(1) move (NOT a ULP drift) — the unambiguous behavior-neutrality violation.

**Lesson.** A "behavior-neutral" claim about zeroing a field is only as good as the SINGLE fission/emission contract it was proven against — when the field is a SHARED source feeding multiple consumers with DIFFERENT contracts (here: same-region χ·νΣf for SN vs sink-region-χ × source-region-νΣf for peierls MG), the inertness must be re-proven for EVERY consumer, not assumed. "Snapshots didn't move" + "no guard ValueErrors" is necessary but NEVER sufficient for behavior-neutrality of a shared data edit: the snapshot suite may exercise only the consumer where the field IS inert, and a guard-error count is blind to silent accuracy regressions in a different solver. Prove behavior-neutrality with a DIRECT old-vs-new value comparison on EACH consuming code path (O(1) move = not neutral), and run the (possibly slow) accuracy-band suites that consume the edited field — do not rely on a fast proxy. → numerical-bug-signatures cross-cutting: H5 (test count is not coverage) + the L20 shared-source hazard (a guard-at-source edit's blast radius includes every consumer of the source, each with its own contract).

**ROOT-CAUSE CORRECTNESS FIX (numerics-investigator, 2026-06-21, branch `feature/field-typed-operator-algebra`, NOT committed).** The QA entry's option (i) — "fix `solve_peierls_mg` to use the SOURCE region's χ" — is the physically-correct resolution and has now been made. The sink-indexed χ was a genuine correctness bug (Mode-2 variable swap: sink index `i` where source index `j` belongs), not just a behavior-neutrality hazard. **Physics ground (structurally independent of the ORPHEUS theory page, which documented the SAME wrong convention — reference contamination, vv §6):** Hébert (2009) *Applied Reactor Physics* Ch.3 Eq. 3.57/3.58 — the fission emission density `Q_fiss(r) = χ_g(r)·Σ_g' νΣ_f,g'(r)·φ(r)` is a single LOCAL quantity at the fission point r; χ shares the source argument with νΣf and φ. The transport kernel `K(r_i,r_j)` is the sole carrier of spatial coupling and never touches χ. **Independent structural witness:** the sibling Variant-α `trajectory_resolvent` solver family already source-indexes χ correctly (`chi_nodes = chi[region_at_node]`, local product `chi_nodes·F_r`), confirming source-indexing is right by a structurally-distinct kernel. **Twin sites fixed (3, all `chi[i]→chi[j]`):** `peierls_nystrom/geometry.py:6416/6426` (`solve_peierls_mg`, all curvilinear: slab-polar/cyl/sph delegate here), `peierls_nystrom/slab.py:313` (native `_build_system_matrices` interior loop), `slab.py:368` (white-BC re-entry). **No re-baseline needed:** the fix is byte-identical (Δ = 0.000e+00) on uniform-χ data — proven by stash-vs-fix on HEAD's equal-χ library; the bug was masked because HEAD gave every region (fissile A + non-fissile B/C/D) IDENTICAL χ (`[1,0]` 2G, `[0.6,0.35,0.05,0]` 4G), so χ_i≡χ_j everywhere. The S10a precursor zeroing B/C/D χ broke that equality and exposed it. **Which test catches it (root-cause level):** the intrinsic-property gate "a non-fissile region's χ must not affect k_eff" — `tests/derivations/` (promoted from `derivations/diagnostics/diag_err063_probe_e_intrinsic_property.py`, `@pytest.mark.l1 @pytest.mark.catches("ERR-063")`, cylinder+sphere negative legs + a positive non-degeneracy leg; Mode-8-safe `np.testing.assert_array_less`). Mutation-verified teeth: reverting `chi[j]→chi[i]` reddens both negative legs (4.7% k_eff move on a non-fissile χ flip) while the positive leg stays green. **Documentation debt (archivist):** `docs/theory/references/peierls_nystrom.rst` Eq. `peierls-mg-operator` (~line 1116) + `tilde B` (~1134) + `solve_peierls_mg` docstring math (geometry.py:6273) all still write the sink-indexed `χ_g(r_i)` and must be corrected to `χ_g(r_j)`; the promoted test's `verifies("peierls-mg-operator")` edge should only be trusted after that correction.

## ERR-064 — SN `compute_keff` omitted the boundary-leakage term: the reported k on a vacuum-bounded eigenproblem was P/A, not the posed problem's eigenvalue — +87% on a leaky slab, +22% on the het vacuum sphere — while the flux converged CORRECTLY, so every flux-based and reflective anchor stayed green

**Date:** 2026-07-03 (#291 symptom of the #259 estimator fragmentation; branch `refactor/k-estimator-unification`; filed at #290 P5 when diffusion's loss-derived denominator made the SN copy's omission visible).
**Module:** `sn` (`orpheus/sn/solver.py::compute_keff`).
**Class:** silent wrong-answer in a REPORTED functional — the iteration converges cleanly to the right eigenvector, then labels it with a non-eigenvalue.

**Failure mode:** #3 missing term (leakage), structurally protected from every existing gate by the estimator-independence of the flux fixed point. `power_iteration` renormalizes the iterate to unit production each outer step, so the k estimate enters only as the `Fφ/k` scaling that renormalization cancels — the flux trajectory and its fixed point are IDENTICAL under any estimator, and `converged()` consumes dk of a (biased but converging) sequence. The wrong k is therefore invisible to: all flux/shape assertions, all reflective k anchors (leakage ≡ 0 ⟹ P/A degenerates correctly), all estimator-shift-invariant equivalence gates (SI≡Krylov, 1-D≡2-D — both legs share the functional), and the particle-balance tests (they recompute P/A and compare TO the reported k — same functional, anti-pattern #8's telescoping cousin). Measured (map-ratio ground truth, noise ≤2e-11): homogeneous 2G vacuum slab w=8 reported **1.8377** vs true **0.9816** (+87.2%, L/A=0.872 — a near-critical leaky slab reported strongly supercritical); het vacuum sphere P0 0.8648 vs 0.7060 (+22.5%); reflective control bias 1.2e-10.

**How it hid.** No vacuum ABSOLUTE-k anchor existed anywhere in tests/sn — the only vacuum-k consumers asserted orderings and banded DIFFERENCES of the reported k (`test_p1_lowers_heterogeneous_keff` band, `test_vacuum_keff_lower_than_reflective` ordering), all of which the biased functional satisfies. The protocol docstring (`numerics/eigenvalue.py`) had stated `k = production/(absorption + leakage)` all along — the implementation violated its own contract, and nothing compared reported-k to an independently posed eigenvalue.

**Fix.** The unified estimator discipline (#259 P1): `k = R_νΣf/(R_Σa + L − E_2n)` with `L` the net outflow through vacuum faces — the face-area integral of the NEW `AngularBoundaryFlux.net_current(face)` (the `Σ_m (Ω·n̂)w_m ψ_m` contraction, single-sourced from the trace space's `omega_dot_n` + partial-current metric), read from the solver-held trace of the last inner solve with a fission-production scale bridge. Reflective faces are a STRUCTURAL zero (skipped, never accumulated) so all-reflective Σ₂=0 problems remain BITWISE the historical P/A — zero re-baseline blast radius outside vacuum problems. vv decision: principled re-baseline (the two P1-effect sphere docstring values re-measured; the band survived — Δ doubled to 2.73e-2, still in (1e-3, 5e-2)).

**Which test catches it.** `tests/sn/eigenvalue/test_keff_estimator_gate.py::TestReportedKeffIsThePosedEigenvalue` — reported k ≡ the converged-fixed-point map ratio `k* = P(Mφ*)/P(φ*)` (structurally independent: consumes only the iteration map + a linear functional, never the A/L/E decomposition), across {vacuum slab, vacuum sphere (pins the 4πR² face-area convention against `volume_measure`), reflective bitwise pin, Σ₂≠0}. Teeth mutation-verified in-file: the leakage-drop mutant re-opens the 87% bias on vacuum AND leaves reflective bitwise green (that asymmetry IS the catcher); the sign-flip mutant is caught by the scale-bridge guard (crash-RED).

**Lesson.** A REPORTED functional whose consumer cancels it (renormalized power iteration) is a class of its own: the solver can be perfectly correct in every state variable while labeling the answer wrongly, and every self-referential or difference/ordering assertion stays green. The committed catcher must compare the functional against the POSED problem's definition (here the iteration map's own Rayleigh quotient at the fixed point) — not against the functional itself in another guise. Corollary of Mode 12: the flux fixed point's estimator-invariance group contains every estimator error; only a gate OUTSIDE that group (functional vs map) sees it.

## ERR-065 — SN/MoC put the (n,2n) emission in the k numerator while every inner solve scales only FISSION by 1/k: the reported ratio equals the posed problem's eigenvalue only when Σ₂=0 or k=1 — −26% at zero leakage on a Σ₂-bearing reflective slab — and MoC's L0 pinned the wrong convention as intended

**Date:** 2026-07-03 (#259 P1 recon finding, ruled R7; sibling of ERR-064 — same "estimator disagrees with its own posed problem" class, different term).
**Module:** `sn` (`compute_keff` numerator `+2Σ₂`) + `moc` (`core.py::compute_keff` same spelling + the L0 `test_n2n_1g_analytical_keff` pinning `(νΣf+2Σ₂)/(ΣC+ΣF+Σ₂)` as the expected value).
**Class:** convention mismatch between estimator and posed problem, plus a self-referential test reference (the L0's "analytical" value was derived FROM the estimator's own convention, not from the eigenproblem's balance — anti-pattern #7/#11 family).

**Failure mode:** #6 convention drift. All four solvers' `solve_fixed_source` treat the (n,2n) emission as a plain (UNSCALED) gain; only fission carries 1/k. The 1G reflective balance is then `(1/k)νΣf + 2Σ₂ + Σs = Σt` ⟹ `k* = νΣf/(ΣC+ΣF−Σ₂)`. The SN/MoC functional `(νΣf+2Σ₂)φ/Σa_xs φ` evaluates to `k* + 2s₂(1−k*)/Σa` — exact only when s₂=0 or k*=1. Measured: reflective slab with the #269 asymmetric-Σ₂ fixture, reported **1.92857** vs true **2.61278** (−26.2%; the closed form `0.78/(0.5185−0.2200)` checks the map-ratio exactly). CP and diffusion were already operator-consistent (CP: `Σt−Σs−2Σ₂` denominator, L1-pinned against a dense eigensolver with `SigS_eff = SigS+2Σ₂`; diffusion: n2n inside the loss action via `IsotropicN2N`).

**How it hid.** Every tight SN k-anchor uses Σ₂=0 library mixtures (the #269 Mode-10 observation — the term contributes exactly 0 on library data), and MoC's lone Σ₂≠0 value test derived its expected k from the SAME convention it was testing — a reference that cannot disagree with the code by construction. The estimator fragmentation (#259) let two conventions coexist because no cross-method law said which was right; the diffusion template's loss-derived denominator (where n2n lands on the removal side automatically) exposed the fork.

**Fix.** R7: operator-consistent everywhere — νΣf-only numerator; the emission `E_2n` enters the denominator NEGATIVELY (`R_Σa + L − E_2n`; `absorption_xs` counts the collision once, the emission is a gain). MoC's L0 expected value re-derived from the posed balance (1.125 → 1.25; principled re-baseline — the old reference encoded the estimator, not the eigenproblem). `compute_production_rate` deliberately KEEPS total physical production (fission + emission) — it is the ERR-052 renormalization scale anchor, a different role from the k numerator (role split documented at both sites). Σ₂=0 arithmetic bitwise unchanged (adding/subtracting exact zeros).

**Which test catches it.** `test_keff_estimator_gate.py::TestReportedKeffIsThePosedEigenvalue::test_reflective_n2n_convention` (reported ≡ map-ratio k* on the Σ₂≠0 fixture; the old convention mutant reds at −26% via `TestGateTeeth::test_old_convention_reds_n2n_leg`) + the re-derived `tests/moc/test_verification.py::test_n2n_1g_analytical_keff` (now posed-problem-derived, k=1.25 exact for any flux shape).

**Lesson.** The k functional is DEFINED by what the inner solve scales by 1/k — estimator and posed problem may not disagree, and "which terms sit in the numerator" is not a style choice. A reference value derived from the code's own convention is self-referential (the L11 structural-independence requirement applies to the DERIVATION of the expected value, not just to solver cross-checks): derive test references from the posed problem's balance. When a term is zero on all library data (Σ₂), its convention errors are Mode-10 invisible — the catcher needs the term genuinely activated (the #269 fixture idiom) AND a reference outside the estimator's own algebra (the map ratio).

## ERR-066 — the SN 1-D adjoint matvec DROPPED the degenerate pure-azimuthal rows: `loss_action_transpose` had no volumetric branch, so `(L+C)ᵀ` on any cylinder quadrature with |μ_x| < eps ordinates (the DEFAULT `product(8,8)` included — 4 | n_phi samples φ = π/2, 3π/2) was silently wrong, while every cylinder gate rode `level_symmetric`, whose point sets have no such ordinate

**Date:** 2026-07-04 (#280 Phase 2.5a — found by code-reading the fold target during the apply-loop unification, instrument-confirmed red before the fix).
**Module:** `sn` (`loss_representation/__init__.py::_OneDimScanWalk.loss_action_transpose`).
**Class:** incomplete transpose — the reverse-mode program omitted an entire branch of its primal (the forward matvec's degenerate volumetric branch), so the coded "transpose" was not the transpose of the coded forward on the degenerate subspace.

**Failure mode:** a missing-branch sibling of #3 (whole rows of Mᵀ absent, not a factor wrong within a row); the blindness mechanism is Mode 7 **in the quadrature dimension** — the nulling ansatz was the RULE CHOICE, not the flux ansatz. The forward writes `m[deg,i] = (denom·ψ − angular_numer)/V` through its degenerate branch (A_downstream = 0, zero faces, no march); the adjoint accumulated neither the `denom·ō/V` diagonal into `psi_bar` nor the `−ō/V` angular-numerator cotangent into `numer_bar`. Measured: G-reciprocity ⟨Aψ,φ⟩_G = 6.923011e+01 vs a mismatched ⟨ψ,Aᵀφ⟩_G — an O(1) relative defect on a 4-cell 2G product-quad cylinder.

**How it hid.** Every cylinder row in the adjoint gate chain — `test_g_adjoint_reciprocity`, the frozen `walk_matvec_cyl_2g` baseline, the removal-form transpose canary — builds on `level_symmetric(4)`, whose direction cosines are bounded away from zero: the degenerate class was EMPTY in every gate, so the missing reverse branch contributed exactly zero everywhere it was tested. Equispaced product rules populate the class whenever 4 | n_phi (φ = 2πk/n_phi hits π/2 iff n_phi/4 is an integer — the default `product(8,8)` qualifies), so production `.H.apply` on a standard product-quad cylinder was reachable-wrong.

**Fix.** The degenerate reversal block in the 2.5a loop-walk frame (a consumer of the shared `_degenerate_positions` trichotomy): slot-local `psi_bar[:, deg, i] += denom·ō/V` and `numer_bar[level][:, within, i] += −ō/V` per degenerate position; `closure.angular_adjoint` reverses the M-M thread's ψ-dependence exactly as for leg ordinates. No face terms (the primal branch has none). Zero participation when the degenerate class is empty — all pre-existing gates bit-identical.

**Which test catches it.** `tests/sn/operators/test_g_adjoint_reciprocity.py::test_g_adjoint_reciprocity_full_block[cyl_product_2g]` (the `_make_cyl_product` builder — red pre-fix at O(1), green post-fix at rel < 1e-12, dev-verified in that order) + the structural trichotomy pin `tests/sn/sweep/core/test_one_dim_loop_walk.py::test_legs_and_degenerates_partition_the_quadrature` (legs ∪ degenerates partition the ordinate set — a two-threshold drift cannot silently re-empty the class).

**Lesson.** A transpose is complete only with respect to EVERY branch of its primal: enumerate the forward's branches (leg march, degenerate volumetric, boundary writeback) and demand a reverse arm per branch — reverse-mode completeness is a per-branch checklist, not a per-function one. Mode 7's "declare what the ansatz nulls" discipline extends to QUADRATURE selection: a gate family whose every row uses one rule family silently nulls the classes that family never populates; the cheap catcher is one builder on the activating rule (here product with 4 | n_phi), not a new oracle. The G-reciprocity gate sees any dropped branch the moment one activating row exists — tolerance tightness was never the issue, coverage of the input lattice was.

## ERR-067 — the SN curvilinear starting-direction (ψ½) block Hilbert metric was `G_sd ≡ 0` (the "ghost metric") — an OPERATOR coefficient mis-installed as the STATE metric, so `A.H = G⁻¹AᵀG` SEVERED the seed→bulk coupling and was a WRONG adjoint for any nonzero seed, while every reciprocity gate fed a present-but-ZERO seed and stayed green

**Date:** 2026-07-06 (the FaceField codim-1 design session; the ghost framing was authored 2026-07-04 in the #280 2.5d carrier work and refuted here — derived by numerics-investigator via the adjoint constraint `AᵀG = GA†` on the dense augmented sphere operator).
**Module:** `numerics` (`spaces/starting_direction_space.py::for_levels` built `np.zeros`) / `sn` (the augmented loss composite's G-adjoint).
**Class:** Mode 12 (invariant-functional gate blindness) — a state's inner-product metric conflated with an angular-integration weight.

**Failure mode.** THREE quantities all vanish at the pole μ=±1 for DIFFERENT reasons, in DIFFERENT objects, and the bug welded them: **M1** the moment/output weight `w_moment = 0` (open Gauss–Legendre has no node at the pole — correctly excludes ψ½ from `φ = Σ wₙψₙ`); **M2** the angular through-flux `(1−μ²)|_pole = 0` (the α-dome endpoint `α_{1/2}=0` — an OPERATOR coefficient in the angular-derivative term, governing the straight-characteristic/triangular structure); **M3** the STATE metric `G_sd` — the inner product, governing the G-adjoint. Installing M2 (=0) as M3 put the seed in `ker G`, so `A.H = G⁻¹AᵀG` (with `G⁺_sd = 0`) zeroed the adjoint seed block `A†_sb = G⁺_sd A_bsᵀ G_b` — the seed→bulk coupling `A_bs ≠ 0` (the Morel–Montry recurrence) went unmatched. Measured production `A.H` reciprocity: **4.6e-16 on a zero seed, 1.3e-2 on a random seed** — a wrong adjoint the instant ψ½ carries data. ψ½ is a first-class radial STATE field (its self-block `A_ss` is a banded radial transport operator `μ∂_r + σ_t`, off-diag ‖71‖ vs the trace's reflection-map ‖2‖), so its Hilbert metric is the bulk's spatial measure restricted to the ray: `G_sd = V_cell`.

**How it hid.** Every G-adjoint reciprocity gate (`test_g_adjoint_reciprocity`, the #282 seed gates, the `full_field_space` metric tests) fed a present-but-**ZERO** ψ½ block — and with a zero seed the coupling term `ψ½ᵀ A_bsᵀ G_b φ_b` is identically zero regardless of `G_sd`, so the seed rows sat inside the reciprocity functional's invariance group: **designed-green at every tolerance, in every regime** (the exact Mode-12 mechanism — annihilation, not sub-floor). The blindness was even DOCUMENTED as correct (a "Mode-12 honest-scope note"), which is how a wrong adjoint shipped behind a green wall. Reachable-wrong the moment the #276/#281 adjoint work routes φ* through `A.H` on a curvilinear (nonzero-seed) problem.

**Fix.** `G_sd = V_cell` (radial cell volume, SPD), installed in `StartingDirectionSpace.for_levels` (mirroring `G_bulk = V_cell·w_n`; the angular factor `w` is a free gauge — a single pole ray has no canonical quadrature weight — gauge-fixed to `V_cell` so the adjoint seed block is the physical backward radial march and bulk/trace observables are bitwise gauge-invariant). The metric is pinned only up to SPD by `AᵀG = GA†` (since `apply_transpose` is the EXACT Euclidean transpose, `‖T−Aᵀ‖ ≈ 3e-16`); `0` is the single forbidden value (the boundary of the SPD cone). **Forward path bit-identical** (the metric is read only by `A.H` + `inner_product`, never the forward sweep/matvec — the #208 trace-metric-install precedent).

**Which test catches it.** `tests/sn/operators/test_starting_direction_metric.py::test_derive_gsd_and_close_mode12` (`@pytest.mark.catches("ERR-067")`) — the dense structurally-independent `G⁺TG` sweep over {0, id, V_cell, V·w} × {zero-seed, random-seed, seed-row-flip}, proving `g_sd=0` is blind AND `V_cell` closes; corroborated on the production path by `tests/sn/sweep/curvilinear/test_282_direct_seed_fixed_point.py::test_mode12_g_reciprocity_catches_a_seed_row_flip`. Mutation-verified teeth: revert `G_sd → 0` reds both; a seed-row sign-flip reds; the corner-gauge `V[-1] → 2·V[-1]` leaves reciprocity + bulk observables green (gauge inert). **Critical gate-design note:** the closure gate carries a **control leg** (unmutated nonzero-seed reciprocity `< 1e-12`) so the red is attributable to the flip — without it, a broken `G_sd=0` baseline (also `> tol` on a random seed) mimics "caught."

**Lesson.** A state DOF's Hilbert metric is set by its OPERATOR ROLE (it must be SPD so the adjoint is non-degenerate), NOT by its output/integration weight — a DOF that is zero-weight in the output functional (`φ`, `k`) can still require a nonzero metric for its adjoint (`output-functional weight ≠ state inner-product metric`). Corollary for Mode-12-closure gates: the blindness is exact precisely on the input the old gate fed (the zero seed), so closing it REQUIRES exercising the previously-nulled input (a nonzero seed) AND a control leg that pins the honest baseline — else the closure test cannot distinguish "the metric now catches it" from "the baseline is broken too." The three-quantities-conflated pattern (integration weight / operator coefficient / state metric all vanishing at one geometric locus) generalises to any augmented block on a measure-degenerate face.

## ERR-068 — `B_a.apply_transpose` output-projected `lawᵀ` onto the outflow rows: the transpose of the row-projected forward `B_face = P_inflow ∘ law` was spelled `P_outflow ∘ lawᵀ`, extracting the vacuum mask law's harmless identity-on-outflow block into a spurious +1 diagonal exactly where the forward is the ZERO map — while every reciprocity fixture rode permutation laws, on which both spellings are bit-identical

**Date:** 2026-07-11 (B.2d d3 — caught by the NEW A2a grid-reciprocity arm on its FIRST run, defect 2.6e-4 on the het-VACUUM coherence sphere; localized by dense discriminator to exactly 4 diagonal trace entries of the (A,A) quadrant — expected `G⁻¹AᵀG` −1 vs `grid.H` −2 = LC's honest −1 PLUS the spurious +1; fixed same session @ `1de9592`).
**Module:** `sn` (`orpheus/sn/operators/boundary.py::SNBoundaryOperator._reflect_trace`, the Euclidean transpose arm).
**Class:** Mode 2 (composition-order swap on a transpose: `(P ∘ law)ᵀ` spelled as `P′ ∘ lawᵀ` instead of `lawᵀ ∘ P`) — hidden by the ERR-063 masking family (a fixture regime on which the wrong and right spellings coincide) and PINNED green by anti-pattern #12.

**Failure mode.** The forward boundary action is row-projected: `B_face = P_inflow ∘ law` — the law maps the full face trace, and only the inflow rows are written. The honest Euclidean transpose is `B_faceᵀ = lawᵀ ∘ P_inflow`: restrict the INPUT to the forward's codomain rows, write the FULL `lawᵀ` image. The implementation instead computed the full `lawᵀ` image and wrote only the OUTFLOW rows (`P_outflow ∘ lawᵀ`). Output projection extracts a law's DIAGONAL block: the realized VACUUM law is the full-face mask (zero-on-inflow ⊕ identity-on-outflow), whose identity block is harmless in the forward (the row projection annihilates it) but leaked under the wrong transpose as a spurious `+1` on every outflow diagonal slot — `B_aᵀ` was `+1` exactly where `B_a` is the ZERO map. Forward solves never call `B_aᵀ`, so no k/flux anchor could ever see it; it was reachable-wrong for every transpose consumer (`grid.H`, reciprocity, the future φ* adjoint chain).

**How it hid.** Off-diagonal permutation laws (reflective, albedo) have NO diagonal block — both spellings select the same entries bit-identically, so every reflective-fixture reciprocity gate was green over the wrong arm (the ERR-063 masking family: neutrality proven on one law-contract is silently false on another). Worse, `test_adjointable_when_all_faces_support` had PINNED the output-projected behaviour as a green snapshot (anti-pattern #12 in its purest form — a pin on an arm no production consumer exercised, certifying the masked regime and actively defending the bug).

**Fix.** (`1de9592`) `_reflect_trace`'s transpose arm masks the input and writes the full image (`masked[sel] = face_in[sel]; out[...] = law.apply_transpose(masked)` — the `(P_sel ∘ law)ᵀ = lawᵀ ∘ P_sel` spelling), with the raw verb licensed by the `adjointable()` TypeGuard + a loud per-face `MissingAdjoint` (spec §39.1; ratchet-clean where the old `getattr` was pyright-invisible). `test_adjointable_when_all_faces_support` rewired onto the defining law (`got = lawᵀ(P_inflow·ψ)`). Forward path untouched BY CONSTRUCTION (the arm has no forward callers) — the wall's kinf gates confirm; `B_b` audited clean (KIND-dispatched, vacuum zero both directions).

**Which test catches it.** `tests/sn/operators/test_psi_half_coupling.py::TestBoundaryUnweld::test_b_a_vacuum_transpose_is_the_honest_zero` (`@pytest.mark.catches("ERR-068")`) — densified vacuum `F ≡ 0` AND `T ≡ 0` on the masking-free regime; born RED against the pre-fix tree, green post-fix, dev-verified in that order. Companions: `test_b_a_transpose_is_the_dense_euclidean_transpose` (the defining law `dense(B_aᵀ) ≡ dense(B_a)ᵀ`, bit-equal on the non-trivial reflective arm) and the end-to-end catcher that found it live — the A2a grid `.H` reciprocity arm in `tests/sn/operators/test_inverse_adjoint_coherence.py` (< 1e-12 + the metric-strip tooth).

**Lesson.** Transpose the COMPOSITION, not the pieces around the wrong side: a row-projected forward `P ∘ law` transposes to `lawᵀ ∘ P` (input restriction, full image written), never `P′ ∘ lawᵀ` (output projection) — and the two coincide exactly on laws with no diagonal block, so a fixture family of permutation laws CANNOT distinguish them. Before crediting a transpose/neutrality gate, enumerate the law/contract classes and check the fixture activates one on which the candidate wrong spellings DIFFER (here one vacuum-fixture row was the whole cost — the third ERR-063-family instance). And a snapshot pin on an arm with zero production consumers certifies the masked regime, not the claim: retire or rewire such pins onto the defining law.

## ERR-069 — the d3 direct-seed carve wired the ψ½ r=R corner datum for vacuum(0) and reflective(B_b) only: a PRESCRIBED inflow's corner went silently ZERO, and the sphere prescribed-inflow solve converged to a wrong fixed point (L2 ≈ 0.21 vs the 2.4e-3 floor) for nine days — because the change was verified exactly on the two regimes where the missing arm is accidentally exact, and the one committed catcher was slow-marked

**Date:** introduced 2026-07-04 (`a29ab2df`, #282 route (a) d3 — the direct starting-direction seed); caught + fixed 2026-07-13 (coupled-block campaign step 7, `1a5e5dc5`) by the tier's FIRST patient slow run since d3.
**Module:** `sn` (`solver._build_fixed_source_rhs`) / `transport` (`RadialCharacteristicField.source_from_angular`).
**Class:** Mode 9 (verified only in the degenerate/coincident regime) compounded by a slow-marked-only catcher — five not-slow walls (6371, 6388, 6385, …) ran green over the wrong fixed point.

**What happened.** Pre-d3, the lagged seed march started its r=R inward
leg from `bc_outer_value = inflow_full[most_inward_global, :]` — the
boundary trace's most-negative-μ ordinate, an O(Δμ) nearest-node proxy
for the seed BVP's Dirichlet datum ψ_in(μ=−1, R) that served ALL three
BC families through one read. d3 retired that read and re-posed the
corner as an explicit given-data slot with per-family arms: vacuum ⇒ 0
(the factory's zero corners), reflective ⇒ the `B_b` corner GAIN
(iterate-side, in the splitting's N). A prescribed inflow is NEITHER —
it is a SOURCE, and no source-side arm existed: the corner stayed the
factory's zero, System B marched from a wrong Dirichlet datum, and the
seed thread fed the wrong ψ½ into every level's Morel–Montry
recurrence.

**Why every gate stayed green.** The §16.D acceptance matrix (sphere
cold-residual, seed-insensitivity, coarse-LS, pure-absorber) was
vacuum/reflective-only — the exact regimes where the missing arm is
zero-exact (vacuum: 0 IS the datum; reflective: B_b supplies it).
Mode-9's signature precisely: the formulation change was FP-invariant
on the verified class and wrong outside it. The one committed catcher
(`test_mms_prescribed_inflow_sphere_activates_redistribution`, whose
band was DESIGNED for this dropped-q.boundary class) is slow-marked, so
the canonical not-slow walls never executed it; the not-slow structural
sibling (`test_sphere_nonvacuum_inflow_honoured_and_redistribution_live`)
checks inflow-honoured + redistribution-live but not the corner datum,
so it stayed green over the wrong converged value.

**Fix.** (`1a5e5dc5`) The inflow-corner law completed to THREE arms,
each datum through its own system channel:
`source_from_angular(..., boundary_trace=)` — the SAME composite rhs's
boundary member — populates each carrying level's inflow corner with
the trace's most-inward `xmax` row (the pre-d3 proxy restored through
the source channel; the half-range Legendre fold to μ=−1 is the named
upgrade seam); `_build_fixed_source_rhs` threads its own boundary
member in. Dormant by DATA off the prescribed class (vacuum/reflective
rhs traces are zero): the 257/0 dormancy battery (byte-compare
walk_matvec / affine-carve / bc_extraction stores, the regression
floor, the r1 reflective-corner row) passed with ZERO re-baselines.

**Which test catches it.**
`tests/sn/operators/test_psi_half_coupling.py::TestPrescribedCornerDatum`
(`@pytest.mark.catches("ERR-069")`) — the NOT-SLOW wiring pin: the
prescribed arm (a random ordinate-asymmetric trace must land on every
carrying level's inflow corner as the most-inward row; the outflow
defect corner stays zero) + the vacuum control (bare-ndarray rhs keeps
zero corners — the dormancy leg). Teeth mutation-verified in-process:
the factory reverted to the d3 state (drop `boundary_trace`) ⟹ RED in
0.14 s; restored ⟹ green. The value-level catcher remains the slow
MMS gate (red→green at the fix, 4/4 file).

**Lesson.** A refactor that re-poses one shared read as per-family
arms MUST enumerate the family — the BC lattice here is
{vacuum, reflective, albedo, prescribed}, and "the two arms my
acceptance matrix exercises" is not an enumeration. Mode-9 discipline
applies to BC-family coverage exactly as to splittings: verify on a
config where the arms DIFFER (a nonzero prescribed trace), not only
where the missing one is accidentally exact. And a bug class whose
only catcher is slow-marked is invisible to every routine wall — when
a slow gate is the designed catcher for a regression class, give the
class a not-slow WIRING pin (the value gate stays slow; the wiring pin
runs everywhere), and run the patient slow tier at every merge gate.

## ERR-070 — the σ_r-fold: realizing the within-group scattering gain as a diagonal σ_r-sweep (`A_wg.solve` with removal σ_t − σ_s0^{g→g} and the isotropic self-scatter gain dropped) treats `Σ_s0·P_iso` as `Σ_s0·I` — exact for isotropic flux, a measured 43% fixed-point error on anisotropic flux, while every isotropic-regime test stays green

**Date:** class identified 2026-06 (#215 investigation — the fold was
never wired into production; the guard docstrings at
`scattering.py:968-971` date from it); formally re-caught and filed
2026-07-26 by the #2 DSA battery's seeded mutation (Phase 3c).
**Module:** `sn` / `transport` (the within-group solve's operator
spelling — the trap any "fold σ_s into the removal" optimization
springs).
**Class:** Mode 9 (exact on the degenerate isotropic-flux regime,
value-wrong outside it) — vv-principles' canonical splitting instance
(test-design mode 9(a)).

**What it is.** The isotropic within-group gain is the rank-structured
`Σ_s0^{g→g}·P_iso` (scatter INTO isotropy from the angular MEAN). The
fold rewrites the within-group solve as a pure-attenuation sweep with
removal `σ_r = σ_t − σ_s0^{g→g}` — algebraically `Σ_s0·I` (scatter
into isotropy from EACH ORDINATE separately). The two coincide iff the
angular flux is isotropic (`P_iso ψ = ψ`); the difference operator
`σ_s0(I − P_iso)` annihilates exactly the isotropic subspace.

**Measured (the seeded reproduction, 3c).** Heterogeneous 2-zone
(A|B) 2g slab, K = 40, S4, vacuum walls, isotropic scattering — flux
anisotropy from the vacuum boundary layers alone: max relative
scalar-flux shift **43.2%** (group 0: 11.4%, group 1: 43.2%) —
4.3e+06× above the D3 FP-identity band (1e-7). #215's own configs
measured 46–56%. On the fully-reflective isotropic box the shift is
IDENTICALLY zero — the designed-green degeneracy that made the fold
look safe.

**Why it never shipped.** The #215 investigation fenced the accessors
(`MaterialXSField.foldable_sigma` / `residual_sig_s`) behind guard
docstrings before any consumer wired them into a sweep; the DSA
low-order build (#2) is their FIRST production consumer — legitimate
there and only there, because the correction→0 property makes the
low-order operator value-safe by construction (a wrong low-order
degrades the rate, never the answer).

**Which tests catch it.**
- `tests/sn/acceleration/test_dsa_rate.py::TestSigmaRFoldCaught`
  (`catches("ERR-070")`) — the seeded fold (realized as folded
  Mixtures: `σ_t − diag Σ_s0`, diagonal zeroed, balance preserved)
  moves the FP 43% ⟹ D3's `ψ_DSA ≡ ψ_SI` band reds on any wiring of
  the fold into the accelerated path; the companion leg proves DSA
  faithfully accelerates even the folded operator (no laundering —
  only the comparison against the TRUE plain FP exposes it).
- `tests/sn/acceleration/test_dsa_rate.py::TestD10RoutingSentinel`
  (`catches("ERR-070")`) — the structural fence: an AST sweep of
  `orpheus/` confirms the foldable accessors' consumers are exactly
  {the definition site, the split layer, the DSA low-order build};
  a new consumer reds at design time, before numerics can be wrong.
  Tooth: a planted `mat_xs.foldable_sigma()` module is flagged.
- `tests/sn/architecture/test_stage_separation.py::
  test_the_sigma_r_fold_is_a_splitting_only_with_its_anisotropic_remainder`
  (`catches("ERR-070")`, added 2026-07-29 @ `9a546640`) — the **algebraic**
  catcher, and the earliest of the three: it constructs the fold's splitting
  `M = (L+C) − Σ_s0·𝕀` and gates the splitting law `A = M − N`, so the defect
  reds at CONSTRUCTION, before any solve runs. Three legs, all measured:
  the honest complement `N = −Σ_s0(𝕀 − P_iso) + B` satisfies the law
  (2.93e-17); `N = 0` reds on an anisotropic state (5.43e-03 relative, 2.63
  absolute); `N = 0` is machine-invisible on an angularly-flat state
  (3.57e-18) — the Mode-9 degeneracy, shipped as a permanent control leg so
  the anisotropic leg's claim to be *what caught it* stays falsifiable. The
  fixture is vacuum-both-faces so `B ψ ≡ 0` (measured 0.0) and the legs
  isolate the projection mechanism alone. A companion test asserts the
  mechanism directly (`Σ_s0·P_iso ≡ Σ_s0·𝕀` on a flat flux, both 2.045) and
  deliberately carries NO marker — it characterises the degeneracy without
  detecting the defect, so a marker there would be a phantom coverage edge.

**Lesson.** A splitting that is exact on a degenerate subspace needs
its FP-invariance gate run on a config that EXITS the subspace
(vacuum/heterogeneous ⟹ anisotropic flux — scattering anisotropy not
required), and the data path that makes the wrong spelling convenient
deserves a structural tripwire alongside the numerical gate: the
fold's 43% is invisible to every isotropic fixture, every balance
check, and every convergence-order study.

## ERR-071 — the composite sweep inverse dropped the rhs's OUTFLOW-trace rows: `(L+C)⁻¹` seeded its boundary buffer from `rhs.boundary`, then let the march clobber the outflow slots — an inverse exact on every physical rhs (outflow rows ≡ 0 there) and SINGULAR on the full composite space, discovered when the P1-DSA GMRES preconditioner stalled at an O(1) true residual with a machine-zero preconditioned residual

**Date:** latent since the composite-carrier sweep (the #208/W-C era);
caught and fixed 2026-07-26 (#2 Phase 3c — the P1-DSA d₁ Krylov
posture excited the kernel deterministically; the end-of-solve
ConvergenceCertificateError made the catch).
**Module:** `sn` (`operators/streaming.py::StreamingCollisionOperator.
_solve_timed_full_field`).
**Class:** Mode 9 in its purest structural form — the operator was
EXACT on the degenerate subspace every physical path lives in (rhs
outflow-trace rows are 0 = 0 defect identities at every builder), and
WRONG on its complement, which only a Krylov iteration's residual
vectors ever visit. No physical fixture could have caught it.

**What happened.** The composite forward `(L+C)` carries the boundary
as a sibling block: inflow rows are identities on the given inflow;
outflow rows are the self-consistency defect `streamed − ψ_out`. The
exact inverse must therefore emit `ψ_out = streamed − rhs_out`. The
sweep seeded its mutable boundary buffer per-face from `rhs.boundary`
(consuming the INFLOW rows correctly — the ERR-069 fix), but the march
then overwrote the outflow slots with the marched outflow, silently
dropping `rhs_out`: `(L+C)⁻¹` mapped every pure-outflow-row rhs to
(essentially) zero.

**How it surfaced.** The DSA-preconditioned GMRES posture (#2) builds
`M = (I + 𝒞)∘(L+C)⁻¹`. With the ℓ≥1 gain active (the P1-DSA arm), the
Krylov residual acquired a pure outflow-row component; measured
`‖M q‖/‖q‖ = 1.07e-15` on that vector — M singular — so full-restart
GMRES stalled at an O(1) TRUE residual while its preconditioned
residual sat at 1e-31, scipy reported not-converged, and the
**end-of-solve certificate** (`_certify_within_group_exit`) refused
the claimed convergence: "the honest equation residual is 1.49" — the
#290-era certificate machinery catching a genuinely new class. The
same singular M existed under the 3b P0 posture and the #200 identity
preconditioner era; nothing excited the kernel until so ≥ 1.

**The diagnosis chain (recorded because the instrument order
mattered).** (1) The corrector was probed directly — linear, healthy;
(2) `(I + 𝒞)` was materialized — smallest |eig| = 1.0, NO kernel;
(3) a per-call spy on the production preconditioner found
`min ‖M q‖/‖q‖ = 1e-15`; (4) capturing THAT q showed a pure
outflow-trace vector with `sweep(q) ≈ 0` — the sweep, not the
corrector. A healthy-component-by-component proof chain that
converged on the one composition seam nobody had gated.

**The fix (the root fix, user-ruled R6 2026-07-27 — four parts).**
(1) *Solve half*: one restore in `_solve_timed_full_field` — after the
march, `boundary_buf[outflow rows] −= seed_boundary[outflow rows]`
(the sign pinned by the round-trip gate — the forward's row is
`streamed − ψ_out`). Bit-inert on every physical path (`−= 0`);
nonsingular M for Krylov. Measured: 2G heterogeneous ℓ≥1 Krylov+DSA
converges 287→12 (vacuum) / 305→16 (reflective). (2) *Transpose
half*: the restore's matrix `E_out` is a diagonal partial identity ⟹
symmetric, so `(Aᵀ)⁻¹ = (S_old)ᵀ − E_out` — the SAME one-site restore
in `solve_transpose`'s output-boundary assembly (the G3 full-composite
reciprocity gates red its absence). (3) *Call-site role conversions*:
every caller that passed an ITERATE trace (stale outflow rows) as
`rhs.boundary` — the `solve_sn` eigenvalue-finalize reconstruction and
the `sweep_once` test helper — now routes through the EXISTING
`AngularBoundarySourceSink.prescribed_inflow` factory (inflow slots
only; outflow rows unrepresentable by construction — the Pattern-4
projection already in the tree, bypassed by raw `from_mesh` casts).
The finalize conflation was invisible to every keff gate (interior
marches off inflow slots only) and caught ONLY by the 2-D reflective
trace-balance gate (defect 8.8e-2 at exact keff — a Mode-12 lesson:
the balance functional sees what the eigenvalue cannot). (4) *Honest
scope*: on a cylinder product quadrature the degenerate pure-azimuthal
rows (`μ_r = 0`, excluded from both selectors) are FREE DOFs of the
composite — the forward is a structural zero row there and the
inverse completes with the identity (seed passthrough); the gate
asserts both halves of that pair explicitly (#284's free-DOF slots,
not a partial-inverse regression). (5) *The SCHEDULED sibling stays
source-subspace-scoped*: `ScheduledInvertibleOperator` (the G-S
`M = (L+C) − B_lower`) interleaves mid-march reflects that consume
the buffer's streamed outflow values BEFORE the end-of-march restore
fires, so its inverse is exact only on `{y : y_out = 0}` (every
production rhs; FP-dust round-trips fine). The full-space completion
needs the restore interleaved per-group with the schedule — deferred
until a full-space consumer exists (a G-S-preconditioned Krylov). No
value-guard (a threshold is arbitrary; exact-zero rejects legitimate
FP-dust `M.apply` images): the honest-scope witness is the W2
off-domain characterization pin (`test_gauss_seidel_reification.py`,
the tripwire that REDS when the completion lands), and the production
catcher for a future off-domain consumer is the same end-of-solve
certificate that caught ERR-071. The W2 fixture itself was a hidden
consumer of the old clobber (`LC.solve(random-full-trace)` was only
"trace-consistent" because the drop erased the rhs's outflow rows) —
fixed to build from an inflow-only rhs.

**Which test catches it.**
`tests/sn/operators/test_sweep_inverse_identity.py`
(`catches("ERR-071")`) — the composite round-trip identity
`(L+C)∘(L+C)⁻¹ ≡ I` on a RANDOM composite with every block populated,
parametrized over {vacuum slab, reflective slab, product-quadrature
cylinder}, plus the pure-outflow-row leg (the exact singular subspace:
the old sweep mapped it to zero), plus the in-process mutation tooth
(an emptied outflow selector reproduces the pre-fix annihilation).
The transpose half:
`tests/sn/operators/test_loss_transpose_solve.py::test_g3_full_field_solve_reciprocity`.
The call-site conflations: the 2-D reflective trace-balance gate
(`test_bc_extraction_2d.py::TestBoundaryResidual2DDrivesToZero`) and
the Q/Σt component-iteration gate (`test_solver_components.py`).

**Lesson.** An inverse operator's contract is the identity on the
WHOLE space, not on the subspace physical data happens to span — and
the round-trip gate `A∘A⁻¹ ≡ I` on a random full-space vector is
cheap, total, and catches every partial-inverse class at once. Wire
it the day the inverse is born, not the day a Krylov method wanders
into the kernel. (And: an exit certificate that re-checks the honest
equation residual converts a silent wrong-answer stall into a loud
refusal — it is the reason this entry exists.) The root-fix corollary:
the drop survived because callers CONFLATED roles — an iterate trace
(state) cast raw into a source slot (given data) — and the correct
projection (`prescribed_inflow`, outflow rows unrepresentable) already
existed; a hand-rolled `from_mesh(trace.values.copy())` beside a
Pattern-4 factory is the same smell as a hand-rolled loop beside an
operator: the cast site is where role errors live, so spell role
conversions through the typed factory, never through raw values.

---

## ERR-072 — `SubgroupOfO3.SO2.is_invariant` certifies a NON-`SO(2)`-invariant rule: the "representative orbit" of a CONTINUOUS group is four rotations `{0°, 90°, 180°, 270°}`, which generate only `C_4`, so every product quadrature with `n_phi ≡ 0 (mod 4)` passes

**Date:** filed 2026-08-02 (pre-carve audit of
`orpheus/numerics/symmetry.py`, ahead of widening `_orbit_closure` to
return its orbit permutation).
**Module:** `numerics` (`symmetry.py:769-779` the generator;
`:665-670` the `SO2` dispatch; `:672-675` the `O2` dispatch, which
inherits it).
**Class:** Mode 12 (invariant-functional / designed-green) in its
group-theoretic spelling — the *checked group* is a proper subgroup of
the *claimed group*, so every measure in the gap is annihilated before
any assertion sees it.

**What it is.** `_orbit_closure` is a **generator-set** check: closure
under every listed matrix implies closure under every product, i.e.
under the group the listed set *generates*. That makes the strategy
sound for `C_n`, `D_nh`, `O_h`, `I_h` — each listed set generates its
group. `_so2_representatives()` returns `[rot_z(kπ/2) for k in
range(4)]`; the first is the **identity** (a vacuous check) and the
group generated is `C_4`, not `SO(2)`. Any `C_4`-invariant rule
therefore certifies as `SO(2)`-invariant.

**Measured.** `Quadrature.product(n_mu=4, n_phi=8)`:
`SO2.is_invariant → True`, while closure under `rot_z(22.5°)` (half
the azimuthal spacing) is `False` and closure under `rot_z(1°)` is
`False` — the rule is exactly `C_8`-invariant. Same `True` for
`n_phi ∈ {4, 8, 12, 16}`; `False` for `n_phi ∈ {2, 3, 5, 6, 7}`, so
the answer is a function of `n_phi mod 4` rather than of the rule's
actual symmetry. The self-inconsistency is visible without any
external reference: the (sound) `C_n` sibling correctly REJECTS
`Cn(16)` on the same measure while `SO2 ⊋ C_16` is accepted —
violating `A ⊆ B` ∧ `B`-invariant ⟹ `A`-invariant.

**How it hid.** Three layers. (1) The docstring pre-authorises the
gap — "*necessary but not sufficient* for general measures, but
**sufficient by construction** for the rules ORPHEUS ships" — which is
FALSE for the product family, one of the four rules shipped. A caveat
that names the risk reads as if the risk had been assessed. (2) The
only `SO2` test in the suite feeds a **1-D** measure, which takes the
`_check_invariance_1d` short-circuit and never reaches the 3-D path:
an instrumented whole-suite run measured `_so2_representatives` called
**0 times across 182 tests**. (3) `is_invariant` returns a bare
`bool`, so a `True` from a `C_4` sample is indistinguishable from a
`True` from the real group — the return type cannot carry the evidence
that would expose it.

**Which test catches it.** None today (the defect is filed, not yet
fixed). The gate that catches it by construction is the
**monotonicity law** — for every asserted lattice edge `A ⊆ B` and
every measure `m`, `B.is_invariant(m) ⟹ A.is_invariant(m)` — which
flagged 68 violations on 11 measures × 19 groups in one loop, of which
48 trace to this entry. The companion is a return-shape change: once
`_orbit_closure` returns its per-operator permutations,
`len(perms) == 4` for a group of infinite order is visibly wrong.

**Lesson.** A finite "representative orbit" certifies **only the group
the sample generates** — compute that group and compare it to the
claimed one, never sample-and-hope. For a CONTINUOUS group the honest
discrete statement is usually a different predicate entirely: a finite
node set on `S²` is `SO(2)`-invariant iff every node lies on the axis
(an off-axis orbit is a circle), so the correct answer is `False` for
every genuine `S²` rule and the `invariance_group=SO2` tags on the
product/GL rules are a claim about the **continuum measure being
discretised**, not about the discrete measure. Two different claims
must not share a predicate name. → `numerical-bug-signatures`
Signature 4 (a quadrature-dependent constant baked into a check) in
its group-theoretic form.

---

## ERR-073 — `_orbit_closure` documents "find a permutation π such that M(nodes)_i = nodes_{π(i)}" but only ever finds SOME same-weight partner per node: the map is never checked for injectivity, so duplicating one node of an `O_h`-invariant rule yields a measure with `M#µ ≠ µ` that certifies invariant

**Date:** filed 2026-08-02 (same pre-carve audit as ERR-072).
**Module:** `numerics` (`orpheus/numerics/symmetry.py:904-954`;
the same defect in `_is_reflection_invariant_1d:619-648`).
**Class:** Mode 12 — the measured functional ("does every image have a
same-weight partner?") has the many-to-one maps inside its invariance
group, so the multiplicity error is annihilated exactly, at every
tolerance, on every input.

**What it is.** Invariance is `M#µ = µ` as **measures**. The loop
computes, for each `i`, some `j` with `‖nodes[j] − M·nodes[i]‖ ≤
100·atol` and `|w_j − w_i| ≤ atol`, and never asserts that `i ↦ j` is
injective (hence, on an equal-size set, bijective). A many-to-one
match satisfies every assertion in the body while certifying a measure
whose pushforward carries different **multiplicities**.

**Measured.** Take `level_symmetric(sn_order=4)` (24 nodes, genuinely
`O_h`-invariant) and append an **exact duplicate** of node 0 (25
nodes, total mass 12.566 → 13.090). `is_invariant(O_h)` returns
`True`. The match map is non-injective for **48 of 48** operators (24
distinct targets from 25 sources). The measure is demonstrably not
invariant: for the first non-identity `M`, `mass at p0 in µ = 1.0472`
vs `mass at p0 in M#µ = 0.5236` — a factor of two. **No tolerance
games are involved**: the duplicate is bit-identical, so no loosening
of `atol` could rescue the check.

**How it hid.** The docstring asserts the strong claim ("find a
permutation π") while the body implements the weak one, and the
`-> bool` return makes the two indistinguishable to every caller and
every test. No test in the suite builds a duplicated-node measure —
the negative cases are a 2-point asymmetric set and a Lebedev-vs-`I_h`
rejection, both of which fail on **partner absence**, never on
**multiplicity**. The bug is also invisible to the four shipped rules,
which are honest permutations: adding the injectivity check leaves
`lebedev(7/17/23)`, `LS(2/4/6/8)` and `product(4,8)` all `True`, so
the fix costs nothing and buys nothing *until* a caller consumes the
permutation — which is exactly what the pending carve does.

**Which test catches it.** None today. The catcher is the duplicate
construction above (`invariant rule + one duplicated node ⟹ must be
rejected`), plus `sorted(perm) == list(range(n))` asserted on the
carve's returned permutation.

**Lesson.** "Every image has a matching partner" is **not** "the map
is a bijection" — the second is what invariance means, the first is
what a nearest-neighbour loop computes, and on a `-> bool` return they
are the same green. Whenever a docstring names a **structure**
(permutation, bijection, isomorphism, partition) that the body only
*implies*, either assert the structure or weaken the docstring; and
prefer **returning the structure** over returning a bit about it,
because a returned permutation makes its own bijectivity assertable
while a `bool` makes it unfalsifiable. (Corollary for orbifold /
stabiliser work: the fixed-point set that makes a node *singular* is
`{i : perm[i] == i}` — derivable from the returned permutation and
underivable from the `bool`, so the structure was load-bearing all
along.)

---

## ERR-074 — `_find_reflections` computed a reflection partner map by bare `argmin` and never checked that the reflection PERMUTES the nodes: on an odd-`n_φ` product rule the stored `σ_x` map was wrong by 0.58 in the direction cosines **and was still an involution**, so the one committed check passed on it, and that map feeds the cylindrical `r = 0` pole continuation

**Date:** filed 2026-08-02 (quadrature-machinery campaign, Q5.0.1 —
found while establishing the fold's preconditions, not by a failing test).

**Where:** `orpheus/numerics/quadrature/directional.py::_find_reflections`
(retired at `a7695148`), consumed via `Quadrature.reflection_index` by
`orpheus/sn/loss_representation/__init__.py` (the `r = 0` pole seed) and
by the specular / reflective BC realizers.

**The bug.** The partner map was

```python
ref[i] = np.argmin( ||R x_i - x_k||²  over k )
```

with **no distance threshold, no injectivity check, and no weight
comparison**. It therefore *asserted* a structure — "these are the
reflection partners" — that it never verified. On a node set the
reflection does not preserve, `argmin` still returns an index: the
nearest node to a point which is **not in the set at all**.

**Why it was reachable.** A product rule's mirror planes sit at
`kπ/n_φ`. `σ_x` is `φ → π − φ`, i.e. a plane at `π/2`, which needs
`k = n_φ/2` — an **integer only for even `n_φ`**. So every odd-`n_φ`
product rule is *not* `σ_x`-closed. `[M]` at `n_φ = 5 / 7 / 9` the
stored axis-0 map was wrong by **0.58 / 0.42 / 0.33** in the direction
cosines. (`σ_y` is the `k = 0` plane and survives at every `n_φ`,
which is why the azimuthal fold is unconditional while the centreline
map is not.)

**Why nothing caught it — the load-bearing part.** The only committed
check was `test_q4_2_reflection_is_an_involution`, asserting
`ref[ref[i]] == i`. `[M]` **the garbage map satisfies it.** A
nearest-neighbour match on a nearly-mirror-symmetric set produces a
self-consistent pairing, so the map is a genuine involution that is
simply not *the reflection*. The test also ran only on even-`n_φ`
products, where the map is correct anyway — so it was doubly blind.

This is the **same failure mode as ERR-042**, one layer down.
`orpheus/geometry/boundary/_specular.py` had already ratified the
lesson for specular BC pairings: it requires **three** independent
checks (involution, measure preservation, inflow→outflow) *because each
passes a table the other two reject*, and its own prose says "the
involution does not [catch it] (it can still be its own inverse)". The
quadrature layer carried **none** of the three.

**Which test catches it.** `tests/numerics/test_quadrature_directional.py`:
`test_q4_5_every_stored_partner_map_really_is_its_reflection` (all three
checks on every stored map), `test_q4_7_odd_n_phi_product_has_NO_x_reflection`
(the axis is absent and the call raises), and
`test_q4_8_certification_drops_an_axis_when_closure_actually_breaks` (a
perturbed node must drop its axis). Teeth proven by monkeypatching back to
the legacy body: **6 red**, with a positive control confirming the mutation
bit first.

⚠ **A gate-design trap this exposed, worth more than the bug.** With only
**even-`n_φ`** rules in its parametrization, the new three-check gate stayed
**GREEN** under the legacy body — every map it inspected happened to be
correct. A gate over "all the rules we ship" is *not* a gate over "the cases
that can break"; the discriminating parameter has to be in the list on
purpose. The odd rows are what give it teeth.

**The fix, and why it cost nothing.** `symmetry._orbit_closure` already
computes the permutation while proving closure, and already requires a
**bijection** with matched positions *and* equal weights (that requirement
is ERR-073's fix). So routing the partner map through it supplies all three
checks for free, and the permutation it returns **is** the partner map. An
axis the measure is not closed under is now omitted, so `reflection_index`
raises rather than returning a wrong answer.

**Lesson.** A map named for a *group action* ("reflection partners",
"rotation index", "octant image") is a claim that some `g ∈ G` **permutes**
the set. Nearest-neighbour matching computes a map; it does not compute that
claim, and the two are indistinguishable on the sets where the claim happens
to hold. Verify the action, or do not name the map after it. And note the
asymmetry that made this survive: an **involution check is cheap and
reassuring and nearly worthless here** — the wrong map was self-inverse. The
check that bites is the one against the *definition* (`R x_i == x_{ref[i]}`),
which is exactly the check nobody writes because it looks tautological.

Cross-refs: **ERR-042** (the same failure for specular BC pairings, one layer
up, already carrying the three-check discipline), **ERR-073** (the bijection
requirement in `_orbit_closure` that this fix consumes), **ERR-072** (the
sibling "a certification claimed and never computed" in the same module).

---

## ERR-075 — the SN realizer returned an AFFINE operator for prescribed inflow and withheld the `BlockRole.BOUNDARY` stamp on the reasoning that the stamp's absence fenced it out of the `B` block: `SNBoundaryOperator._face_laws` never filtered on the stamp, so the source was delivered through `−B` as well as through `q_∂` — a 2× inflow on source iteration and a hard `ConvergenceCertificateError` on Krylov, both invisible because the law is not a registered `BC` kind and no test ran a full solve with a DECLARED law

**Date:** filed 2026-08-05 (affine-boundary-source campaign, P3 — found by
asking `B.apply(0)` while auditing a *different* claim, not by a failing test).

**Where:** `orpheus/sn/boundary/realizer.py` (the `PrescribedInflow` arm) and
`orpheus/sn/boundary/angular.py::IncomingSourceOperator` (both retired at P3);
the un-filtering consumer is
`orpheus/sn/operators/boundary.py::SNBoundaryOperator._face_laws`.

**The bug.** The affine boundary law is `γ₋ψ = L γ₊ψ + q`. The realizer tier's
contract is to realize `L`. For prescribed inflow `L = 0` and the whole content
is `q` — but the arm returned an operator whose `apply` **ignored its input and
returned `q`**:

```python
def apply(self, psi_out):
    return self.source.evaluate((self.n_inflow,) + psi_out.shape[1:])
```

That is an affine map declared as a `LinearOperator`. `[M]` at `q = 2.5`:
`A(0) = 2.5` and `A(2x) − 2A(x) = 2.5`.

**The load-bearing error is not the affine operator — it is the FENCE ARGUMENT.**
Three docstrings (`numerics/operator.py`, `geometry/boundary/_realizer.py`,
`geometry/boundary/_bound_compat.py`) stated as the design of record that the
rank-0 affine source "is deliberately NOT stamped `BlockRole.BOUNDARY` — it is
the boundary *source* `q.boundary`, not a linear boundary operator `B`", and the
campaign plan drew the conclusion that it therefore could not reach `B`. It
could. `_face_laws` is:

```python
return {face: self.sn_mesh.bc[face]
        for face in self.sn_mesh.angular_trace.layout.faces}
```

— every face, no `block_role` test anywhere. `SNBoundaryOperator` carries the
stamp *itself* and applies all of its face laws, so an unstamped leaf is not
excluded from anything. **The stamp is a classification consumers may query; it
was never a gate.**

**Measured consequences** (heterogeneous 2-group scattering slab, GL-8, identical
composite bulk on every leg; `γ₋(xmin)` read off the converged answer, 12 dp):

| configuration | SI | Krylov |
|---|---|---|
| declared inflow, both channels live | `5.000000000000` | ⛔ raises |
| declared inflow, source channel disabled | `2.500000000000` | ⛔ raises |
| vacuum mesh + hand-supplied `q_∂` | `2.500000000000` | `2.500000000000` |
| vacuum control | `0.000000000000` | `0.000000000000` |

`|φ_double − φ_vac| / |φ_single − φ_vac| = 2.000000` exactly. On Krylov it is
worse and older: an affine `A(x) = A_lin(x) − c` breaks the Arnoldi relation
`A V_k = V_{k+1} H_k`, so GMRES's tracked residual is meaningless and
`_certify_within_group_exit` raises `‖Aψ − q‖/‖q‖ = 1.718`. That path had been
UNUSABLE with a declared prescribed inflow since long before the source channel
existed.

**Why it survived.** One real fence, misread as two. `prescribed_inflow` is not a
registered `BC` kind (`sn/mesh/augmented_mesh.py`, #189), so no production driver
can install the law and only tests construct it — and **no test ran a full solve
with a declared law.** The gates that existed stopped at
`_build_fixed_source_rhs`: they correctly verified the RHS receives `q` and were
structurally blind to `B` also delivering it. This is `verification-user-path`
(`docs/theory/verification/principles.rst`) failing one commit after that
doctrine landed — the gates travelled part of the user's path and stopped short
of the solve.

**The catcher.** `tests/sn/solve/test_declared_inflow_reaches_the_rhs.py::
test_the_declared_boundary_law_holds_on_the_answer[source_iteration|krylov]`
(`catches("ERR-075")`) — asserts `γ₋ψ|_f == q_f` on the converged answer. The
boundary condition is a DEFINITION, so the gate needs no reference solver and no
discretization assumption, and it separates the three outcomes exactly
(`5.0` / `2.5` / `0.0`). Exactness is path-dependent and measured: SI copies the
source into the slot (bit-exact); Krylov reaches the trace through its iterate
(18 ULP at 2.5, gated by `assert_array_almost_equal_nulp(nulp=64)` — a ULP budget
rather than an `rtol`, so it cannot quietly widen into the 1×/2× gap).

**⭐ The generalizable lesson — a "it is fenced, so it cannot happen" argument is
a claim about a CONSUMER and must be measured at the consumer.** Both fences here
were read off the PRODUCER: the realizer returns an unstamped operator, the `BC`
registry has no key. One of those (the registry) happened to be a genuine
consumer-side fact. The other was an *inference* about what `SNBoundaryOperator`
would do with an unstamped leaf, and the code it was an inference about is
fourteen lines long and was never opened. The check costs one line —
`B.apply(zero)` — and a linear operator's answer is `0` by definition, so the
probe needs no oracle. **Whenever a design note says "X is excluded because it
lacks Y", grep for who reads Y.** If nobody does, X is not excluded; and the
absent-`Y` note will read, to every later reader, exactly like a safety argument.

**⭐ Mutation-classification corollary (anti-pattern #18, measured here).** The
obvious mutation — reinstate the affine operator — reddened 17 gates, but it
breaks LINEARITY, so most of those reds see the broken law rather than the
delivery count and must not be credited. The in-class mutations are `q → 2q` in
the source channel (6 reds) and `L := IdentityOperator` (10 reds), both perfectly
linear. The second is the discriminator: `B(0) == 0` still HOLDS under it, so a
leaf-linearity gate is blind and only the converged-trace gate bites. A
superposition gate ("φ is affine in `q`") is a *provable* non-catcher — a doubled
delivery is `q → 2q`, which is still exactly affine in `q` (Mode 12).

Cross-refs: **ERR-069** (the sibling prescribed-inflow delivery bug — the d3
corner datum went silently zero, and also survived because it was verified only
on the regimes where the missing arm is accidentally exact), **ERR-068** (the same
`B_a` block, a transpose-side projection error that rode fixtures where two
spellings coincide), **ERR-047** (the `q ∈ Γ₋` contract this campaign's source
channel now carries), **ERR-063** (a "this change is inert" claim that held for
the contract it was proven against and failed for a second consumer).
