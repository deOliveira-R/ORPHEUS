# Numerics Investigator — Memory Index

Slim index. Behavioral lessons live in `lessons.md` (read FIRST each
session). This index holds only (2) live in-flight state and (3) durable
diagnostic references. Bug-campaign play-by-play lives in the topic
`*.md` files this index points to — do NOT inline it here.

## 1. Lessons (read first)

- [lessons.md](lessons.md) — 17 diagnostic-cascade lessons (L1–L17).
  The spine: never guess, isolate in cascade order, and a single mesh /
  flat flux / homogeneous case / two-probe agreement proves nothing.
  L1 cascade order · L2 curvilinear redistribution is the prime suspect ·
  L3 two-quadrature signed-error gate for rank-N closures · L4
  convergence-rate fingerprints (err·n^p table) · L5 paper-floor vs
  code-bug discriminator · L6 curvilinear matvec needs a NON-FLAT
  per-ordinate hand reference · L7 a direction-dependent per-ordinate
  moment must reach the GLOBAL frame before the angular reduction · L8
  the project's own theory page can be the contaminated reference · L9 a
  "hang" may be a fixture cost — bound the solver first · L10 "diverges
  with refinement" + a discarded library info-flag = an unconverged inner
  solve, not a discretization bug · L11 measure the residual r=Aψ−q for a
  ρ-honest stop, not the increment ‖Δψ‖ · L12 an offline-isolated error
  is THE floor only after an end-to-end swap + silent control + amplify ·
  L13 "diverges/crashes for LD (spectator trailing axis)" = a greedy
  `(Ellipsis, *idx)` fancy-index mis-targeting axes (fix: pin leading axes
  `(slice,slice,*idx)`; non-LD tests are blind) · L14 curvilinear `(L+C).solve`
  is a direct inverse for slab+CYLINDER, seed-LAGGED for SPHERE (MMS-blind) ·
  L15 to rule PRINCIPLED-vs-REGRESSION on a seed/angular-closure re-pose sweep
  ANGULAR N at fixed fine mesh (NOT h at fixed N — two closures share the (h→0,
  N→∞) limit but differ at fixed N); + a carve that GROWS a Krylov composite
  must resize `restart` from the composite `to_flat`, not the bulk formula ·
  L16 to compare two ANGULAR QUADRATURES differing in pole/seed treatment,
  build a standalone scheme-faithful driver gated to production, reference =
  fine-N + a cross-quadrature CONTAMINATION GUARD (MMS is seed/quadrature-
  blind), validate the new pole handling with the per-ordinate flat-flux
  residual · L17 a STATE DOF's Hilbert metric is NOT its angular-integration
  weight; when apply_transpose is the EXACT Euclidean transpose (T=Aᵀ,
  measured) a block metric is GAUGE-FREE (any SPD) — `AᵀG=GA†` only forbids
  DEGENERACY; G_block=0 is worse than blind (breaks the adjoint for
  nonzero-block data); closing a Mode-12 blindness needs the metric AND
  nonzero-block gate probes.

## 2. Active / in-flight state

**None.** Every SN / curvilinear / Peierls campaign this agent diagnosed
is **merged to main** (git-verified: the campaign commits are ancestors
of HEAD) — the ERR-026 family (#195/#196/#168/#193/#229/#9), the
matvec/sweep unify (#206), the eager-registry "hang" (#212), the
pole-cell NaN (#209), the LD diffusion limit (ERR-061 #240), the Peierls
χ source-index (ERR-063 #257 S10a), the affine flux typing (#208), and
the S9 LD boundary-slope verdict. Behavioral lessons are in `lessons.md`;
the V&V record is in the `vv-principles` `error_catalog.md` + the SN
theory page.

> Merge-status in memory goes STALE. ALWAYS reconcile a "resume X" /
> "OPEN #NN" note against `git merge-base --is-ancestor <hash> HEAD` and
> `gh issue view <NN>` before acting — never trust a frozen "NOT landed".

**Open investigations with a diagnostic-of-record** (Peierls/CP rank-N
family — no active work; pick up only if asked): #123 (two-quadrature
rank-N gate), #128 (QMC angular quadrature), #132 / #100 (rank-N Marshak
mode-0/mode-n normalization), #129 (planar-limit physics), #170 (cylinder
bare-critical). Their breadcrumb files are in section 3.

## 3. Durable reference (reusable diagnostic recipes)

- [coupled_block_operator_numerics.md](coupled_block_operator_numerics.md)
  — NUMERICS ASSESSMENT of coupled 2×2 blocks `[[A_AA,A_AB],[A_BA,A_BB]]` for the planned
  `CoupledOperator` (SN ψ½ + DSA + multiphysics; branch `refactor/sn-walk-unification`,
  pre-DSA). The ψ½ system has TWO partitions: ray/bulk WITHIN `(L+C)` is block-TRIANGULAR
  `A_BA=0` (measured `A_sb=A_st=0` EXACT) → DIRECT forward substitution (#284); the full
  `A=(L+C)−S−B=M−N` is the OUTER block-GS, `A_BA` lives in the LAGGED `S` (`S_sb=0.183`;
  `S_bs=0` = ψ½ zero moment weight), `ρ(M⁻¹N)=0.371<c=0.4`. **WRAP-first** (welded sweep =
  exact inverse 3.5e-16; extract-to-dense = principled 5.5e-16, NOT bit-identical; naive
  dense M⁻¹ on arbitrary rhs = O(1) garbage — must preserve the sweep's inflow/seed
  row-contract). Folded ρ=0 strictly beats lagged (edge-extrap #282 diverges ρ≈70; dense LU
  is the oracle, not a rebuilt lagged path). The dagger-biproduct substrate ALREADY EXISTS
  (`FullFieldSpace`+`BlockRole`+`OperatorSum.H`+`GreenOperator` splitting); solve strategy =
  a FAMILY keyed by SPECTRAL CHARACTER (triangular-direct / splitting-iterative /
  preconditioned-Krylov / dense-direct-ELLIPTIC — diffusion's `MatrixInverseOperator`, since
  the Neumann series diverges for fine-mesh elliptic); diffusion `(J⁺,J⁻)` = the Schur
  precedent; 3 block-maths stay distinct (`BorderedOperator` REJECTED, saddle-point deferred
  #294) ⟹ CoupledOperator = a VIEW/strategy over the biproduct, NOT a monolith. ρ-honest stop
  (`evaluate_residual` exists, unused as metric) is a coupled requirement (L11). Diagnostics
  `diag_coupled_0{1,2}_*.py` (6 green). Backs L14.
- [starting_direction_metric_gauge_derivation.md](starting_direction_metric_gauge_derivation.md)
  — DERIVATION-of-record + recipe: the SN curvilinear starting-direction (ψ½)
  block Hilbert metric `G_sd`. Verdict: `apply_transpose` is the EXACT Euclidean
  transpose (T=Aᵀ, 3.6e-16) ⟹ `A.H=G⁺AᵀG` honest for ANY invertible G ⟹ **G_sd
  is GAUGE-FREE** (`AᵀG=GA†` only forbids DEGENERACY; any SPD works). The ghost
  `G_sd=0` is the ONE forbidden value — worse than blind: BREAKS the adjoint for
  nonzero seed (1.3e-2, unmatched `A_bs` coupling). Recommend **V_cell** (radial
  volume, matches `G_bulk=V·w_n`; angular w = sole gauge d.o.f.). The `(1−µ²)=0`
  ghost argument confuses the angular-INTEGRATION weight with the STATE metric.
  Dense unit-vector-probe recipe (check `T==Aᵀ` vs numpy transpose = independent
  ground). Backs L17. Diagnostics `derivations/diagnostics/diag_gsd_0{1,2,3}_*.py`
  (17 green); forward stays bit-identical under the install (#208 trace precedent).
- [radial_characteristic_carrier_level_position_key.md](radial_characteristic_carrier_level_position_key.md)
  — VERDICT (gates the RayOp/A_BB carve): the ψ½ direct-march `cells_view(p_idx)`
  is NOT a p_idx-vs-level bug. `RadialCharacteristicSpace.levels` are LEVEL
  POSITION indices (`enumerate(raw)`), the SAME coord that keys slots + is the
  sweep p_idx + indexes `level_ordinates_list` — and `seed_levels` (gate) ≡
  `space.levels` (key validator) share ONE `radial_characteristic_levels` source,
  so cells_view never crashes AND hits the right slot for ANY carry config. The
  loop's `level` (walk mu_level_idx; None for sphere) is the WRONG key — the
  "fix" pass-level would crash the sphere. Empirically sphere→(0,) always,
  cyl→() never carries (R12a τ_raw∈{0,1}). Carve keys by p_idx. Diagnostic
  `diag_p_idx_vs_level_radial_characteristic.py` (34 green, -O-safe, KEEPER).
- [glob_vs_gl_spherical_quadrature_study.md](glob_vs_gl_spherical_quadrature_study.md)
  — DESIGN-study verdict + recipe: Gauss-Lobatto vs Gauss-Legendre angular
  quadrature for spherical SN (would GLob's nodes-at-μ=±1 turn the M-M seed
  ψ½ into an ordinary weighted ordinate?). Verdict: GLob tracks GL at a
  bounded **~1.2× error penalty** at resolved N (N≥8, N>L), regime/c/aniso-
  insensitive — AFFORDABLE. The μ=−1 GLob node → **τ_0=0** → recurrence
  singular → MUST be straight-char (Carlson march), so no clean production
  point-swap. Method = standalone scheme-faithful driver gated to production
  (1e-11) + fine-N reference w/ cross-quadrature contamination guard (MMS is
  blind). Backs L16. Driver `scratch/experimental/glob_sphere_study/`,
  diagnostics `derivations/diagnostics/diag_glob_0{1..5}_*.py` (33 green).
- [curvilinear_inverse_seed_taxonomy.md](curvilinear_inverse_seed_taxonomy.md)
  — is 1-D curvilinear `(L+C).solve` an honest `SweepOperator`? SLAB yes; **CYLINDER
  QUADRATURE-dependent** (level-symmetric DEAD-seed trivially yes / **product SEED-LAGGED,
  cold err 0.575** — corrected 2026-07-05; the "α-dome telescopes" mechanism was a
  mis-attribution); SPHERE was seed-lagged, FIXED by route (a) `a29ab2d`. **#280 Phase 2.5b
  ruling: the product-cyl lag is FEASIBLE to retire by a PURE-DIAGONAL fold** (κ=dA_w[m0]·
  c_in[m0] into the m0 cell diagonal; POC single-pass = M⁻¹ 5e-16; fixed point machine-identical
  (M⁻¹ to 5e-16) so keff/MMS/matvec gates don't move; change site `_run` 4091–4101, matvec unchanged; LS
  inert, degenerates downstream, #229 clamp defines κ not an obstruction). Mode-7 MMS-blind.
  Backs L14. Diagnostics `derivations/diagnostics/diag_280_cyl_{product_seed_lag,fold_feasibility}.py`.
- [cp_matrix_density_and_sphere_conservation.md](cp_matrix_density_and_sphere_conservation.md)
  — CP [P] is structurally DENSE (0 exact zeros in every geom/tau config →
  #226 name is storage-agnostic `MatrixInverseOperator`; sparse interface-current
  variant = OPEN #56, not built). Probe: isolate [P] via
  `compute_pinf_group(sig_t_g)`, vacuum→P_cell / white→P_inf, reciprocity =
  `diag(SigT·V)@P` symmetric, conservation = rowsums. **Bonus OPEN DEFECT**:
  spherical CP violates the proven row-sum=1 above tau~3 (rowsum saturates at
  1/pi; grows with refinement, flat in n_quad_y), invisible to the tau≤2 suite.
- [issue_208_flux_displacement_residual_typing_debug_value.md](issue_208_flux_displacement_residual_typing_debug_value.md)
  — the convergence-diagnostic method catalogue. FluxDisplacement (Δψ
  tangent — `contraction_ratio` >1 diverges / ≈1 stalls / <1 healthy,
  `where_largest`, sign-oscillation) answers "WHERE is convergence
  lagging"; AngularResidual (r=Aψ−q rate-density — `balance_map`,
  `boundary_vs_interior_split`, `relative_to`, `as_dsa_source`) answers
  "WHERE is the equation unsatisfied". (Backs L11.)
- [curvilinear_tau_clamp_vs_pole_floor.md](curvilinear_tau_clamp_vs_pole_floor.md)
  — characterization of record: the sphere pole-CELL spatial closure is
  O(h) at r→0 (mostly an MMS midpoint-vs-shell-volume-avg artifact +
  literature-accepted inherent first-order); INVISIBLE to the
  volume-weighted L2 gate → WONTFIX-closure + an L∞/pole gate. The
  diag_01..31 probe ladder. (Backs L12.)
- [sn_space_angle_discretization_coupling.md](sn_space_angle_discretization_coupling.md)
  — verdict "space & angle discretization are coupled but needn't be":
  Cartesian fully separable (cross-term≈0); curvilinear genuinely coupled
  via the M-M angular thread injecting into the spatial scan source. The
  separability-probe recipe (diag_sep_*).
- Open rank-N / Peierls breadcrumbs (diagnostics-of-record for the OPEN
  issues in §2): [direction_n_quadrature_baseline.md](direction_n_quadrature_baseline.md)
  (#123) · [frame_5_qmc_quadrature.md](frame_5_qmc_quadrature.md) (#128) ·
  [issue_100_class_b_mr_mg.md](issue_100_class_b_mr_mg.md) +
  [issue_132_augmented_nystrom.md](issue_132_augmented_nystrom.md) +
  [issue_132_cylinder_hebert.md](issue_132_cylinder_hebert.md) (#132/#100) ·
  [issue_129_planar_limit.md](issue_129_planar_limit.md) (#129) ·
  [peierls_greens_variant_alpha_decision.md](peierls_greens_variant_alpha_decision.md)
  + [peierls_greens_phase1_closeout.md](peierls_greens_phase1_closeout.md)
  (the Variant-α reference-solver chain).
