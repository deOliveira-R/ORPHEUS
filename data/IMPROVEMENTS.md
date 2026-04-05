# Data Package — Improvement Tracker

Central registry of improvements for cross-section data infrastructure.

## Tracking Number Format

`DA-YYYYMMDD-NNN` where DA = Data, YYYYMMDD = session date, NNN = sequence.

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context

---

## DONE — Implemented and Documented

### DA-20260405-001 — Mixture utility properties

**Status**: DONE  
**Commit**: 443c790  
**Sphinx**: docs/api/data.rst (autodoc)

Added properties to `Mixture` dataclass:
- `absorption_xs` — fission + capture + (n,alpha) + (n,2n) out
- `total_scattering_xs` — P0 row sum (in + out)
- `in_scattering_xs` — in-group elastic (P0 diagonal)
- `out_scattering_xs` — out-of-group (P0 off-diagonal sum)

Eliminates 6+ copies of the same absorption formula across solver modules.

### DA-20260405-002 — Shared per-cell XS assembly (CellXS)

**Status**: DONE  
**Commit**: e79a363  
**Sphinx**: docs/api/data.rst (autodoc)

Created `data/macro_xs/cell_xs.py` with:
- `CellXS` dataclass (sig_t, sig_a, sig_p, chi per cell)
- `assemble_cell_xs(materials, mat_ids)` function

Replaces identical per-cell XS extraction loops duplicated in CP, DO,
SN-1D, and MOC solvers.

---

## OPEN — Not Yet Implemented

### DA-20260405-003 — Production-weighted fission spectrum

**Priority**: Low | **Effort**: Small

`compute_macro_xs` uses the first fissile isotope's chi spectrum.
Should use production-weighted average:
`chi_mix = sum_i(nu_i * Sigma_f_i * chi_i) / sum_i(nu_i * Sigma_f_i)`.
Currently acceptable for single-fissile problems (UO2) but incorrect
for MOX or mixed assemblies.

### DA-20260405-004 — h2o_properties viscosity TODO untracked

**Priority**: Low | **Effort**: Small

`data/materials/h2o_properties.py:38` has an untracked TODO about using
IAPWS viscosity for all calls instead of pyXSteam.  Related to the NaN
fix via IAPWS fallback documented in memory.  Should decide: either
switch entirely to IAPWS or keep pyXSteam as primary with IAPWS fallback.

### DA-20260405-005 — Scattering matrix Legendre order consistency

**Priority**: Low | **Effort**: Small

`Mixture.SigS` stores P0, P1, P2 Legendre orders but most solvers only
use P0.  The DO 2D solver uses up to PL.  No validation that the
requested order is available in the Mixture.
