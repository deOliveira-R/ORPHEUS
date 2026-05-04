# Orphan Equation Classification

Classification grid for the 78 orphan equations under the four
resolution classes:

- **A**: Should-have-test (write a `@pytest.mark.verifies(...)` test)
- **B**: Documented-only (definitional / derivation-step / asymptotic)
- **C**: Stale (remove)
- **D**: Already tested elsewhere (add label to existing test)

## Trajectory-resolvent (`docs/theory/trajectory_resolvent.rst`)

| Label                                                   | Class | Resolution                                                                                                                                                  |
| ------------------------------------------------------- | :---: | ----------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `peierls-greens-defining-bvp`                           |   B   | Defining BVP statement — definitional setup, no isolated test                                                                                               |
| `peierls-greens-A1-split`                               |   B   | Algebra step: split `t = t_bar + t_h` — derivation step                                                                                                     |
| `peierls-greens-A5-specular`                            |   B   | Sanchez (A6) homogeneous-BC kernel form — definitional                                                                                                     |
| `peierls-greens-L0`                                     |   B   | Trajectory-geometry definition `L_0 = sqrt(R^2-r^2(1-mu^2))-r·mu`                                                                                            |
| `peierls-greens-Lp`                                     |   B   | Trajectory-geometry definition `L_p = 2 R mu_surf`                                                                                                          |
| `peierls-greens-mu-surf`                                |   B   | Trajectory-geometry definition `mu_surf(r,mu)`                                                                                                              |
| `peierls-greens-trajectory-integral`                    |   B   | First-leg trajectory integral definition `F(r_i,mu_q)`                                                                                                       |
| `peierls-greens-bounce-period-integral`                 |   B   | Bounce-period integral definition `B(mu_surf)`                                                                                                              |
| `peierls-greens-surface-fixed-point`                    |   D   | Verified by `test_v_alpha1_surface_fixed_point_solves_to_q_over_sigma_t`                                                                                    |
| `peierls-greens-T-mu-surf`                              |   B   | Closed-form `T(mu_surf) = 1/(1-exp(-Σt L_p))` — geometric series sum (definitional)                                                                          |
| `peierls-greens-function-architecture`                  |   D   | The `psi = F + e^{-Σt L_0} T B` closure is the architecture; verified by every Variant α numerical test (`test_peierls_greens_function_solver.py`)          |
| `peierls-greens-V-alpha-1`                              |   D   | Verified by `test_v_alpha1_*` (4 tests) in `test_peierls_greens_function_symbolic.py`                                                                       |
| `peierls-greens-k-inf`                                  |   B   | k_eff = k_inf identity — already exercised by `test_v_alpha1_overall_pass` (V_α1 implies this) — but k_inf def is a textbook identity, mark documented      |
| `peierls-greens-T00-integrand`                          |   D   | Verified by `test_v_alpha2_T00_matches_hebert_via_matrix_path`                                                                                              |
| `peierls-greens-V-alpha-2`                              |   D   | Verified by `test_v_alpha2_*` (4 tests)                                                                                                                     |
| `peierls-greens-V-alpha-3`                              |   D   | Verified by `test_v_alpha3_g_h_vanishes_at_alpha_zero`                                                                                                       |
| `peierls-greens-bounce-sum-alpha`                       |   B   | Bounce-sum recurrence step — derivation algebra                                                                                                              |
| `peierls-greens-T-alpha`                                |   B   | Closed-form `psi_surf = α B/(1 - α e^{-Σt L_p})` — geometric series sum                                                                                     |
| `peierls-greens-mg-source`                              |   B   | Multi-group source definition (textbook identity)                                                                                                            |
| `peierls-greens-mr-trajectory-segments`                 |   B   | `r(s)^2 = r_i^2 - 2 r_i mu_q s + s^2` — analytic-geometry identity                                                                                          |
| `peierls-greens-mr-piecewise-tau`                       |   B   | Piecewise-tau definition (multi-region)                                                                                                                     |
| `peierls-greens-garcia-convention`                      |   B   | Convention conversion (notation bridge, factor of ½)                                                                                                         |
| `peierls-greens-fixed-source-iteration`                 |   B   | Fixed-source iteration scheme (algorithm description)                                                                                                        |
| `peierls-greens-unification-resolvent`                  |   B   | `T = (I - S)^{-1}` cross-domain frame statement                                                                                                              |
| `peierls-greens-slab-bounce-period`                     |   B   | `L_period^slab = 2L/|mu|` geometry identity                                                                                                                  |
| `peierls-greens-slab-asym-monodromy`                    |   B   | Slab asym monodromy matrix `S(α_L, α_R, τ)` definition                                                                                                       |
| `peierls-greens-hollow-sph-outer-only-resolvent`        |   B   | Outer-only resolvent rank-1 form (definitional building block)                                                                                               |
| `peierls-greens-sanchez-A6`                             |   B   | Sanchez 1986 (A6) kernel transcription (literature reference)                                                                                                |

## Peierls main page (`docs/theory/peierls.rst`)

| Label                          | Class | Resolution                                                                                                                  |
| ------------------------------ | :---: | --------------------------------------------------------------------------------------------------------------------------- |
| `peierls-boltzmann`            |   B   | Boltzmann transport equation — pure-definitional (already tagged `boltzmann` is documented; this is its peierls-page twin) |
| `peierls-integral-form`        |   B   | Integral form of the transport equation — definitional                                                                       |
| `peierls-3d`                   |   B   | 3D point-kernel form — definitional reduction                                                                                 |
| `peierls-slab-foundations`     |   B   | Slab Peierls reduction — definitional reduction step                                                                          |
| `peierls-cyl-foundations`      |   B   | Cylinder Bickley-Naylor reduction — definitional reduction step                                                              |
| `peierls-sph-ps1982-foundations` |   B | PS-1982 sphere reduction — definitional reduction step (the reduction itself is class-B; PS-1982 cross-check tests the result, not this label) |
| `peierls-bc-general`           |   B   | Sanchez 1986 (A3.a) general BC parametrisation — definitional                                                                 |

## Peierls Nyström (`docs/theory/peierls_nystrom.rst`)

| Label                                  | Class | Resolution                                                                                          |
| -------------------------------------- | :---: | --------------------------------------------------------------------------------------------------- |
| `peierls-slab-polar`                   |   B   | Slab polar-form recap — derivation step                                                              |
| `peierls-mg-operator`                  |   B   | Multi-group operator form recap — derivation step                                                    |
| `peierls-vacuum-bc-flux`               |   D   | Verified by `test_row_sum_matches_analytical_uniform_source` (slab/sphere classes) — add label      |
| `peierls-vacuum-bc-row-sum-gate`       |   D   | The row-sum gate is exactly what `test_row_sum_matches_analytical_uniform_source` enforces           |
| `peierls-vacuum-bc-slab`               |   D   | Slab analytical reference; verified at row-sum and elementwise — add label to `TestSlabKernelRowSum`|
| `peierls-vacuum-bc-cylinder`           |   D   | Cylinder vacuum analytical reference; verified by `TestCylinderKernelRowSum` etc.                   |
| `peierls-vacuum-bc-sphere`             |   D   | Sphere vacuum analytical reference; verified by `TestSphereKernelRowSum` etc.                        |
| `peierls-white-bc-slab`                |   D   | Wigner-Seitz identity 1/Σ_t for slab; verified by `test_flux_equals_1_over_sigma_t`                 |
| `peierls-rank-n-jacobian-derivation`   |   B   | `J_n^+` definition — derivation step                                                                |
| `peierls-rank-n-P-esc-moment`          |   B   | `P_esc^(n)` moment definition — derivation step                                                      |
| `peierls-half-range-inner-products`    |   B   | `<f,g>_L = ∫f g dμ`, `<f,g>_M = ∫f g μ dμ` — definitional inner products                            |
| `peierls-change-of-basis`              |   B   | `M_nm = <ψ_n^M, φ_m^L>_M` change-of-basis definition — derivation step                              |
| `peierls-M-rank-1`                     |   B   | `M^(1) = √2/2` numerical evaluation — already covered by definition of `peierls-change-of-basis`     |
| `peierls-M-rank-2`                     |   B   | `M^(2) = ((√2/2, √6/6), (0, √3/3))` numerical evaluation — derivation step                           |
| `peierls-WM-WL-asymmetric`             |   B   | `W_M = B^μ W_L` matrix-form derivation step                                                          |
| `peierls-class-b-Jn-canonical`         |   B   | Canonical `J_n^+(r_i)` form — equivalent to `peierls-rank-n-jacobian-derivation`                    |
| `peierls-class-b-hebert-closure`       |   B   | Hébert (1−P_ss)^{-1} closure formula — derivation step                                              |
| `peierls-class-b-Pss-homogeneous`      |   D   | `P_ss(Σt,R) = (1-(1+2τ_R)e^{-2τ_R})/(2τ_R^2)` — verified by `test_v_alpha2_Pss_matches_hebert_via_polar_path` |
| `peierls-cyl-Pss-derivation`           |   B   | Cylinder `P_ss^cyl` integral derivation — derivation step                                            |
| `peierls-specular-R-formula`           |   B   | `R_spec = (1/2) M^{-1}` closed-form — derivation step                                                |
| `peierls-specular-M-tridiagonal`       |   B   | Tridiagonal matrix entries — derivation step (numerical evaluation)                                  |
| `peierls-specular-multibounce-formula` |   B   | `K_bc^{spec,mb} = G·R·(I-T R)^{-1}·P` operator form — derivation step                                |
| `peierls-specular-T-matrix`            |   B   | Sphere T_mn closed-form — derivation step                                                            |
| `peierls-bc-operator`                  |   B   | `S_bc(r) = (T_bc q)(r)` operator definition — derivation step                                        |
| `peierls-svd`                          |   B   | SVD of K — definitional (math identity, no test)                                                     |
| `peierls-operator-form`                |   B   | `Σ_t φ = T_vol q + S_bc` operator form — derivation step                                             |
| `peierls-operator-factorisation`       |   B   | `T_bc = G_∞ ∘ R_∞ ∘ P_∞` factorisation — derivation step                                             |
| `peierls-factored-kernel`              |   B   | `K_bc = G·R·P` factored form — derivation step                                                       |
| `peierls-tensor-G-definition`          |   B   | Escape tensor G^i_n definition (proportionality) — derivation step                                  |
| `peierls-tensor-P-definition`          |   B   | Response tensor P^n_j definition (proportionality) — derivation step                                |
| `c-in-remapping`                       |   B   | `c_in(μ_emit)` Jacobian remapping definition — derivation step                                       |
| `mode-conservation-target`             |   B   | Naive per-mode conservation target — counterexample/derivation step                                  |

## Collision-probability (`docs/theory/collision_probability.rst`)

| Label                            | Class | Resolution                                                                                                |
| -------------------------------- | :---: | --------------------------------------------------------------------------------------------------------- |
| `e1-decomposition`               |   D   | E_1(z) = -ln z - γ + R(z) — verified by `test_E1_log_decomposition_notation_bridge` and Atkinson tests   |
| `peierls-sphere-equation`        |   B   | Sphere observer-centred polar form — derivation step                                                      |
| `peierls-sphere-nystrom`         |   B   | Sphere Nyström assembly — derivation step                                                                  |
| `peierls-sphere-ray-optical-depth` |   B   | Ray optical-depth definition — definitional                                                                |
| `peierls-sphere-G-bc`            |   B   | Surface-to-volume Green's function definition — derivation step                                            |

## Singular-eigenfunction (`docs/theory/singular_eigenfunction.rst`)

| Label                | Class | Resolution                                                                                                                                |
| -------------------- | :---: | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `case-method-eq40`   |   D   | X-function (Atalay Eq 40) — covered by `test_case_method_x_function.py`                                                                  |
| `case-method-eq42`   |   D   | z_0 (Atalay Eq 42) — covered by `test_atalay_z0_table1_isotropic` (already verifies `atalay-eq42-extrapolated-endpoint` — eq42 is the same)|
| `case-method-eq46`   |   D   | Slab criticality eq 46 — covered by `test_atalay_table2_slab_*` (slab) and `atalay-eq46-slab-eq54-sphere-parity` already exists           |
| `case-method-eq5`    |   B   | Validity bound `c ≤ 1 + 1/(3 f_1)` — definitional precondition                                                                            |
| `case-method-eq54`   |   D   | Sphere criticality eq 54 — covered by `test_atalay_eq54_sphere_vacuum_isotropic` and parity-flip tests                                    |

## Discrete-ordinates (`docs/theory/discrete_ordinates.rst`)

| Label        | Class | Resolution                                                          |
| ------------ | :---: | ------------------------------------------------------------------- |
| `vacuum-bc`  |   B   | Vacuum-BC operational definition `ψ_n^in = 0`. The semantics are exercised by every vacuum-BC test (`test_vacuum_keff_lower_than_reflective` etc.); this label is the one-line operational rule, not a verifiable identity. Mark documented (definitional). |

## Summary

- A: 0
- B: 60
- C: 0
- D: 18

The dominance of class-B reflects the structure of the trajectory-
resolvent and peierls-Nyström pages: they are derivation narratives.
Each `:label:` is a numbered step. The verifiable claims are the
**capstone identities** (V_α1, V_α2, V_α3, vacuum-BC closed-form
references, P_ss closed form, e1-decomposition) which are class **D**
— they have rich existing test coverage that needs the label added.

The classification respects:

- Cardinal Rule 1 (correctness): no test fabricated to remove an
  orphan; every class-D label maps to an existing test that
  genuinely verifies the labelled claim.
- `vv-principles` Anti-pattern 5: the V_α SymPy tests already verify
  the converged-to value (algebraic identity), not just convergence
  rate.
- `vv-principles` Anti-pattern 7: each class-D test uses
  structurally-independent evidence (closed-form reference for
  vacuum-BC; SymPy two-path integration for V_α2).
- `algebra-of-record`: the V_α tests are the canonical "Branch 1
  State 1A" verification — `derive_*()` returns dict, foundation-
  tagged tests pin them.
