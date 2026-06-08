---
name: phase5a-moment-consuming-scatter-derisk
description: Methodology — a bit-exact "two paths agree" claim is SOUND only with a structurally-independent ground; the ORPHEUS Y_0^0=1 (unnormalized harmonics) load-bearing fact lets moment-consuming scatter equal the full-angular path with no rescale. Phase 5a LANDED.
metadata:
  type: project
---

# Phase 5a moment-windowing de-risk — the reusable methodology

**Status:** Phase 5a "angular-windowing" carve LANDED on origin/main (commits
`93807aa`→`63719a2`): the 2-D SI iterate's bulk moved from full per-ordinate
`AngularFlux (N,ng,nx,ny)` to reduced SH moments `HarmonicMomentField`. De-risked
BEFORE any production edit; verdict was bit-exact at the S level. This note keeps
the two durable lessons; the step-by-step recipe is now in the landed code +
`docs/theory` (`:ref:sn-angular-windowing`).

## LOAD-BEARING FACT — ORPHEUS uses UNNORMALIZED real harmonics (Y_0^0 = 1.0)
`Y[:,0,0]=1.0` (`orpheus/numerics/basis/spherical_harmonic_basis.py`), NOT
`1/√(4π)`. Therefore `M(ψ)[0,0] = Σ_n w_n·1·ψ_n` is the SAME einsum as
`AngularFlux.integrate_angular` (`einsum("n,ngij->gij", w, ψ)`). Consequences any
moment-windowing / moment-consuming carve depends on:
- `Phi[0,0] == integrate_angular` BIT-exactly; `k = Y_0^0 = 1.0` exactly.
- The scalar flux reads off `Phi[0,0]` with **NO `/k` rescale**. (A normalized-
  harmonic codebase would need `√(4π)`; ORPHEUS does not. Forgetting/adding the
  rescale is the convention-drift bug class this fact pre-empts.)
- Producer-side `/W` (Pattern 7): isotropic+n2n source assembled then `/sum_w` at
  the apply boundary; aniso reconstructed via
  `HarmonicMomentReconstruction.from_Y(S.Y).apply(out_moments)/sum_w`.

## METHODOLOGY LESSON — bit-exact agreement is SOUND only with a structural ground
The de-risk reached ULP=0 on Q2/Q3 (`S_apply_from_moments(M(ψ)) ==
S.apply(TFF(ψ))`) GENUINELY — because the carve is literally `Σ_ℓ R_ℓ Λ_ℓ (M(ψ)[ℓ])`
vs `Σ_ℓ R_ℓ Λ_ℓ M_ℓ(ψ)` and `M_ℓ(ψ)` IS the ℓ-block of `M(ψ)`: same einsum, same
operands, same reduction order. **But ULP=0 alone proves only self-consistency — a
shared Y-normalization bug in kernel+candidate would have passed Q2/Q3 silently.**
The claim is sound ONLY because Q2b added a STRUCTURALLY-INDEPENDENT ground:
candidate vs from-scratch Bell & Glasstone hand reconstruction (raw einsums,
`(2ℓ+1)=3` and `/W` written by hand, NOT copied from the kernel) → rel 3.4e-16
(~1.5 ULP, different reduction order). That is what proves the production kernel is
mathematically correct, not just internally consistent.

**General rule (vv-principles §1, L11):** when two paths agree bit-exactly because
they share the einsum/operands, that is the WEAKEST possible cross-check (procedural,
not structural independence). For ANY bit-exact "two paths agree" de-risk: add one
probe whose math is hand-written from the textbook (different reduction order,
different intermediate naming) and assert ~few-ULP, NOT ULP=0. Bit-exact-where-
shared + few-ULP-where-independent is the honest signature; bit-exact everywhere is
a tautology waiting to hide a shared upstream bug (cf. ERR-032: two antiderivatives
agreeing at 1e-39 via a shared wrong identity).

**Evidence (if re-running):** `derivations/diagnostics/diag_p5a_moment_consuming_scatter.py`
— 6 tests, config `_cartesian_2d_p1_aniso_het_si` LS-S4 2G/4G P1 het fuel|mod.
Q4 non-degeneracy guard (aniso l≥1 max 0.495, P0 scatter ASYMMETRIC
|sig0−sig0.T|≈0.1–0.18) makes Q3 an actual ERR-002-sensitive cross-check, not a
flat-flux no-op.
