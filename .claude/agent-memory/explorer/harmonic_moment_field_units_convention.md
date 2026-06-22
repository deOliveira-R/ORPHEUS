---
name: harmonic-moment-field-units-convention
description: Why HarmonicMomentField.UNITS == ScalarFlux UNITS (not angular-flux) — the no-prefactor SH convention (Y_0^0=1, weights sum to 4π) means the sr of angular flux cancels in projection, so a stored moment carries scalar-flux units. Durable physics context; bug history ERR-039/ERR-051.
metadata:
  type: project
---

Durable physics-convention context for the SN moment machinery. The
file:line specifics this once carried are LANDED (the flagged-WRONG docstring
is now corrected; `HarmonicMomentField.UNITS == SCALAR_FLUX_UNITS`) — re-derive
exact lines via Nexus/grep if needed. The CONVENTION below is the durable part.

**The SH convention (single source of truth = `SphericalHarmonicBasis`):**
- No-`4π/(2ℓ+1)`-prefactor real spherical harmonics. `Y_0^0 = 1` exactly.
  Addition theorem `Σ_m Y_ℓ^m Y_ℓ^m = P_ℓ(Ω·Ω')`.
- Quadrature weights SUM TO 4π (sphere rules). So weights carry `sr`; `4π` in
  code is `4π steradians` (an integral over the unit sphere), NOT a pure
  number — decisive for the units argument.
- Forward projection `M`: `φ_ℓ^m = Σ_n w_n Y_ℓ^m(Ω_n) ψ_n` (the `Σ_n w_n`
  integrates over angle → the `sr` of angular flux cancels against the `sr` in
  `w_n = dΩ`).
- Reconstruction `R`: `q_n = Σ_ℓ (2ℓ+1) Σ_m Y_ℓ^m φ_ℓ^m`. The `(2ℓ+1)` is the
  addition-theorem factor, NOT a Hilbert adjoint: `R ≠ M*` — they differ by the
  moment-space Gram `g_C = diag(4π/(2ℓ+1))` (`R = 4π·g_C⁻¹·S₀`, `M* = g_C·S₀`).
  The four operators `{S₀, M=w·S₀, M*=g_C·S₀, R=(2ℓ+1)·S₀}` are SEPARATELY TYPED
  so conflation is structurally unrepresentable — that typing is what fixed the
  bug history below.

**φ_0^0 = scalar flux EXACTLY:** since `Y_0^0=1`, `φ_0^0 = Σ_n w_n ψ_n = φ`
(no 1/4π, no 1/Y_0^0). `HarmonicMomentField.scalar_flux()` returns the
`(0,0)` moment directly and agrees bit-exactly with
`AngularFlux.integrate_angular()`.

**UNITS verdict:** every stored `φ_ℓ^m` has the SAME units for all ℓ, equal to
ScalarFlux's `1/(cm²·s)` (per group bin). The `sr` cancels in the angular
integral; the `(2ℓ+1)` is dimensionless and lives on R; the moment-space metric
`4π/(2ℓ+1)` lives on the SPACE, not the stored value. So
`HarmonicMomentField.UNITS == ScalarFlux.UNITS` (both are angle-integrated flux).

**Bug history (verified accurate):** ERR-039 + ERR-051 — the `(2ℓ+1)`/`4π`
factor was inconsistently applied; a Wave-0 "ugly fix" removed `(2ℓ+1)` from
`apply_transpose` but installed the WRONG operator (bare `S₀` mislabeled as the
Hilbert adjoint). The fix was to give the four operators (S₀, M, M*, R) separate
types living in ONE place (`SphericalHarmonicBasis` + `SphericalHarmonicSpace`)
so the conflation cannot recur. `HarmonicMomentField` is the first typed-Field
consumer of that space; it is the natural windowed-SI iterate type (the moment
tensor the scattering `M` already produces).
