# Homogeneous Infinite-Medium Solver — Improvement Tracker

Central registry of ALL bugs, improvements, and features for the
homogeneous solver.  **This is the single source of truth.**

## Tracking Number Format

`HO-YYYYMMDD-NNN` where:
- `HO` = Homogeneous module
- `YYYYMMDD` = session date when the item was created
- `NNN` = sequential number within the session
- `00000000` = item predates the tracking system

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context
- **WONT**: decided against, with rationale

## Where TODOs Live

TODOs exist in exactly TWO places:
1. **This file** — every item with its tracking number, status, and context
2. **In code** — at the exact location where the fix goes, with the
   matching tracking number (e.g., `# TODO HO-20260403-001: ...`)

No other files should contain TODOs.  If you find one, it must be
consolidated here and given a tracking number.

---

## DONE — Implemented and Documented in Sphinx

### HO-20260403-001 — Sphinx theory chapter for homogeneous model

Documented in `docs/theory/homogeneous.rst` (1146 lines, zero warnings).
Covers: Boltzmann derivation, multi-group discretisation, matrix form,
(n,2n) treatment, scattering convention, 1G/2G/4G analytical solutions
with full worked derivations, cross-section pipeline (sigma-zero,
interpolation, macroscopic summation), spectrum physics (1/E,
Maxwellian, fast), power iteration algorithm, convergence properties,
flux normalisation, two example problems with auto-generated plots,
verification table, solver comparison.

### HO-20260403-002 — API reference page

Documented in `docs/api/homogeneous.rst` (automodule directive).

---

## OPEN — Not Yet Implemented

### HO-20260403-003 — Fission spectrum averaging for mixed fissile isotopes

**Priority**: Low | **Effort**: Small
**Code location**: `data/macro_xs/mixture.py:147`

Currently, the fission spectrum `chi` is taken from the **first fissile
isotope** in the mixture (simplification documented in Sphinx).  For
mixtures with multiple fissile isotopes (e.g., MOX fuel with both
Pu-239 and U-235), the fission spectrum should be the
production-weighted average:

```
chi_mix = sum_i (nu_i * Sigma_f,i * chi_i) / sum_i (nu_i * Sigma_f,i)
```

This is not urgent because ORPHEUS currently uses single-fissile-isotope
fuels (UO2 with U-235 dominant), but it would be needed for MOX or
thorium fuel cycles.

### HO-20260403-004 — Temperature-dependent comparison plot

**Priority**: Low | **Effort**: Small
**Code location**: `docs/theory/homogeneous.rst` (new plot directive)

Add a plot showing how the neutron spectrum changes with moderator
temperature (e.g., 294K vs 600K vs 900K).  This would illustrate:
- Doppler broadening of the thermal peak
- Shift of thermal peak energy
- Hardening of the spectrum at higher temperatures

Useful for pedagogical purposes and connects to the Doppler reactivity
coefficient.

### HO-20260403-005 — Sensitivity of k_inf to boron concentration

**Priority**: Low | **Effort**: Small
**Code location**: `docs/theory/homogeneous.rst` (new subsection)

Add a study showing k_inf vs boron concentration (0–4000 ppm) for
the PWR-like mixture.  This would demonstrate:
- Linear reactivity worth of boron
- Effect on thermal flux depression
- Connection to chemical shim in PWR operation

### HO-20260403-006 — SymPy derivation for two-group explicit formula

**Priority**: Medium | **Effort**: Small
**Code location**: `derivations/homogeneous.py`

The current `derive_2g()` uses SymPy's generic `charpoly` and `solve`
to find eigenvalues.  The Sphinx documentation now shows the full
intermediate derivation (A^-1, M=A^-1*F, quadratic formula).  A
dedicated derivation that produces the intermediate LaTeX steps (not
just the final k_inf) would back the new documentation and could be
included in the generated RST fragments.
