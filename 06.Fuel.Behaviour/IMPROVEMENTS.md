# Fuel Behaviour — Improvement Tracker

Central registry of ALL bugs, improvements, and features for Module 06 (Fuel Behaviour).

## Tracking Number Format

`FB-YYYYMMDD-NNN` where FB = Fuel Behaviour, YYYYMMDD = session date, NNN = sequence.

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context

---

## IMPL — Implemented, Sphinx Pending

### FB-20260401-001 — DAE-to-ODE restructuring (algebraic stress solver)

**Status:** IMPL  
**Date:** 2026-04-01  
**Files:** `fuel_behaviour.py`, function `_solve_stress()`

**Problem:** MATLAB uses `ode15s` with a mass matrix to solve a mixed ODE/AE system
where stress components are algebraic variables (zeros on mass matrix diagonal).
`scipy.integrate.solve_ivp` does not support mass matrices.

**Solution:** Removed stresses from the state vector entirely. Instead,
`_solve_stress()` builds and solves the `3*(nf+nc)` linear system of stress
equilibrium + strain compatibility + boundary conditions at every RHS call.
The ODE state vector contains only temperatures, fission density, swelling,
and creep/plastic strains (190 variables vs MATLAB's ~280).

**Validation:** Open-gap phase matches MATLAB to 0.006°C at t=1 day. Gas pressure
matches to 0.04%. Gap closure detected at 2.85 years.

### FB-20260401-002 — Closed-gap BC4 rewrite (displacement-based constraint)

**Status:** IMPL  
**Date:** 2026-04-01  
**Files:** `fuel_behaviour.py`, function `_solve_stress()`, BC4 section

**Problem:** Original closed-gap BC4 (hoop strain compatibility across gap) used
MATLAB's differential form `d(eps_h)/dr = (eps_r - eps_h)/r`, which requires
dividing by the gap width. When using the algebraic approach with undeformed
geometry, the fabrication gap width (100 μm) was used instead of the roughness
(6 μm), causing a 17× error in the stress gradient and 10× error in contact
pressure.

**Solution:** Replaced with a displacement-based gap-width constraint:
`clad_r0_in * (1 + eps_h_clad(1)) - fuel_r0_out * (1 + eps_h_fuel(nf)) = roughness`.
This formulation:
- Is physically transparent (deformed gap width = roughness)
- Is linear in the stresses
- Avoids the `1/gap_dr` amplification
- Eliminates dependency on `gap_depsh` (strain offset at closure)

**Result:** Contact pressure 40.7 MPa vs MATLAB 39.8 MPa (2.2%). Gas pressure 6.92
vs 6.93 MPa (0.1%).

### FB-20260401-003 — Gap closure event with safe fallback

**Status:** IMPL  
**Date:** 2026-04-01  
**Files:** `fuel_behaviour.py`, function `_gap_closure_event()`

Gap closure event uses one-shot cache from the last RHS call for the deformed gap
width. Falls back to computing from scratch with try/except if cache is stale.
Terminal event with `direction=-1` (triggers only when gap closes).

---

## OPEN — Not Yet Implemented

### FB-20260401-004 — Closed-gap stress transient settling

**Status:** OPEN  
**Date:** 2026-04-01  
**Priority:** Low

Shortly after gap closure (~3 yr), the Python transitions from negative to positive
clad hoop stresses faster than MATLAB (~0.2 yr vs ~2 yr). This is inherent to the
algebraic-stress-within-ODE approach: the MATLAB DAE solver naturally handles the
stiff stress-creep feedback, while the Python's operator-split approach resolves it
more abruptly. Both converge to the same steady-state stresses at EOC.

Not worth fixing unless sub-year accuracy during the closure transient is needed.

### FB-20260401-005 — Sphinx theory documentation

**Status:** OPEN  
**Date:** 2026-04-01  
**Priority:** Medium

No `docs/theory/fuel_behaviour.rst` exists. Should document:
- The 1D radial thermo-mechanical model (heat + stress + swelling + creep + plasticity)
- The DAE-to-ODE restructuring approach
- The algebraic stress solver (linear system construction)
- Open-gap vs closed-gap boundary conditions
- Gap closure event detection
- Validation against MATLAB at t=1 day and EOC

### FB-20260401-006 — Derivation scripts for stress solver

**Status:** OPEN  
**Date:** 2026-04-01  
**Priority:** Medium

No SymPy derivation scripts exist in `derivations/` for the fuel behaviour module.
Key equations needing backing:
- Cylindrical stress equilibrium and compatibility (→ linear system coefficients)
- Hooke's law compliance coefficients `C_self = 1/E`, `C_cross = -nu/E`
- Von Mises stress from deviatoric components
- Prandtl-Reuss flow rule for plastic strain rate
