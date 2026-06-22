---
name: issue-208-flux-displacement-residual-typing-debug-value
description: Reusable convergence-diagnostic signatures from #208 affine flux typing — ‖Δψ‖ false-converges as c→1 (1/(1−ρ) understatement); residual r=Aψ−q is ρ-honest. Method catalogues for FluxDisplacement (Δψ tangent) vs AngularResidual (rate-density). #208 LANDED.
metadata:
  type: project
---

# #208 affine flux typing — the durable convergence-diagnostic signatures

**Status:** #208 (subsumes #201) LANDED on origin/main (commits `8c2f355`→
`9504146`, 2026-06-08): `FluxDisplacement` (affine difference space V) + flux-add
gate (`flux+flux→TypeError`, `flux−flux→FluxDisplacement` kept) + SI
`contraction_ratios` diagnostics + box-7 typed `from_balance` residual eval. This
note keeps only the **reusable numerical signatures** the design exploration
surfaced — the recommendation/decision scaffolding is now history (see the commits).

## SIGNATURE 1 — `‖Δψ‖` false-converges as c→1 (the SI stall-masking trap)
The SI iterate increment `Δψ = ψ⁽ⁱ⁾−ψ⁽ⁱ⁻¹⁾` and the equation residual `r = Aψ−q`
carry DIFFERENT convergence information; related by one preconditioned step
`Δψ = L⁻¹·r̃`. They are NOT interchangeable as a stopping test:
- **`‖Δψ‖` (flux units) UNDERSTATES the true error by `1/(1−ρ)`.** SI contraction
  `ρ ≈ max Σ_s/Σ_t`; `‖Δψ⁽ⁱ⁾‖ ≈ ρ·‖Δψ⁽ⁱ⁻¹⁾‖`; true error `e⁽ⁱ⁾ = Δψ/(1−ρ)`
  (geometric tail). At c=0.99, ρ≈0.99 → "converged at tol" is actually ~100·tol
  from the fixed point. This is the canonical stall-masking-as-convergence trap
  (Adams-Larsen 2002). **Fingerprint:** a solver that reports converged but whose
  answer is wrong by a c-dependent factor; the factor GROWS as c→1.
- **`‖Aψ−q‖` (rate-density) is ρ-HONEST** — measures distance from the solved
  equation, independent of the iteration scheme; does NOT shrink artificially as
  ρ→1. Cost: +1 matvec/iter (or reuse Krylov's). Needs `‖r‖/‖q‖` (relative) to be
  tolerance-portable (its magnitude scales with Σ_t·V).
- **Cross-family:** SAME amplification mechanism as
  [[kinf_homogeneous_curv_mg_inner_tol_amplification]] (ρ/(1−ρ) × inner_tol) and
  [[krylov_restart_truncation_bug]] (info-flag discard). When a "converged" answer
  is wrong, suspect ρ-blind stopping before suspecting the operator.

## SIGNATURE 2 — diagnostic method catalogue by typed role (which question each answers)
Two distinct typed roles; same `(N,ng,nx,ny)` shape + quadrature metric, DIFFERENT
algebra. Use as a menu when building convergence/balance diagnostics.

### FluxDisplacement (Δψ tangent vector — "WHERE is convergence lagging?")
- `contraction_ratio(prev)` → ratio≈ρ; **>1 DIVERGES** (wrong fixed point / unstable
  scheme — the ERR-026 / Krylov-restart family); **≈1 STALLED** (c→1 slow mode,
  curvilinear/reflective); **<1 healthy**. THE single highest-value method (turns the
  ρ-blind `‖Δψ‖` honest). LANDED as SI `contraction_ratios`.
- `true_error_estimate = ‖Δψ‖/(1−ρ)` → catches Signature-1 false-convergence.
- `where_largest` / per-cell-group-ordinate convergence MAP → localised
  non-convergence: pole-cell resonance ([[ordinate_scan_pole_cell_resonance]] NaN
  family), material-interface slow modes, one group lagging.
- sign-oscillation detector → iteration-scheme bug (BC-ordering, G-S reflect-order
  ERR-056) vs discrete-operator bug (the sign-mixed-oscillating fingerprint,
  numerical-bug-signatures §2-D fingerprint).

### AngularResidual (r = Aψ−q rate-density — "WHERE is the equation NOT satisfied?")
- `balance_map` (per-cell/per-group balance violation) → **the typed form of the
  Signature-1 per-ordinate flat-flux residual probe.** Exposes per-ordinate balance
  failure that conservation (telescoping) HIDES (vv-principles H3): e.g. the
  curvilinear pole 22%-error SI gives while global balance holds
  ([[issue_196_phase_g_step2_diagnostic]]).
- `boundary_vs_interior_split` (FREE from the typed `TimedFullField(bulk, boundary)`
  composite) → discriminates BC-realiser bugs from interior-streaming bugs — the
  diagnostic the open [[vacuum_bc_eigenvalue_divergence]] investigation needed.
- `relative_to(source) = ‖r‖/‖q‖` → tolerance-portable honest stopping test.
- `reduction_ratio(prev)` → ρ-INDEPENDENT convergence health (honest as c→1, unlike
  FluxDisplacement's contraction ratio) — distinguishes "converging slowly" from
  "stalled at wrong fixed point".
- `as_dsa_source` → the genuine forward consumer (#2 consistent DSA): the typed
  residual IS the DSA low-order-solve substrate.

## Code seams (may have moved post-#208; verify before relying)
- SI conv test was `numerics/iteration.py:~484` `‖psi−psi_prev‖/max(norm,1e-30)` —
  FLAT ravelled Euclidean (`np.linalg.norm(to_flat())`), IGNORES the angular
  quadrature metric (unlike `Field.l2` which uses `space.inner_product`). The
  flux-add gate kept `‖Δψ‖` as the SI stopping test (bit-identity constraint) and
  ADDED the typed residual as a diagnostic — did NOT switch the criterion.
