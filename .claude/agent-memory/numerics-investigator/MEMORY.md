# Numerics Investigator — Memory Index

One line per entry. Behavioral lessons live in `lessons.md` (read FIRST each session);
campaign detail lives in the topic files this index points to — never inline it here.

## 1. Lessons (read first)

- [lessons.md](lessons.md) — L1–L18, the diagnostic-cascade spine: never guess, isolate in
  cascade order, and a single mesh / flat flux / homogeneous case / two-probe agreement
  proves nothing. Index of the headlines:
  - L1 run the cascade in order · L2 curvilinear redistribution is the prime suspect
  - L3 two-quadrature signed-error gate for rank-N closures · L4 convergence-rate fingerprints
  - L5 paper-floor vs code-bug discriminator · L6 curvilinear matvec needs a NON-FLAT reference
  - L7 a per-ordinate moment must reach the GLOBAL frame before the angular reduction
  - L8 the project's own theory page can be the contaminated reference
  - L9 a "hang" may be fixture cost — bound the solver first
  - L10 diverges-with-refinement + a discarded info-flag = an unconverged inner solve
  - L11 measure the residual r=Aψ−q for a ρ-honest stop, not the increment ‖Δψ‖
  - L12 an offline-isolated error is THE floor only after swap + silent control + AMPLIFY
  - L13 a greedy `(Ellipsis,*idx)` fancy-index mis-targets axes under a spectator axis
  - L14 curvilinear `(L+C).solve` seed-lag is QUADRATURE-dependent (MMS-blind)
  - L15 sweep ANGULAR N at fixed mesh to rule principled-vs-regression on a closure re-pose;
    a carve growing a Krylov composite must resize `restart` from the composite `to_flat`
  - L16 comparing two angular quadratures needs a standalone scheme-faithful driver +
    fine-N reference with a cross-quadrature contamination guard
  - L17 a STATE DOF's Hilbert metric is NOT its angular weight; T=Aᵀ ⟹ the block metric is
    GAUGE-FREE (any SPD); `G_block=0` is worse than blind
  - L18 to adjudicate a LABELING/ORDERING degeneracy the instrument is the operator's own
    SYMMETRY GROUP — MMS is exactly blind when its ansatz uses only the class invariants

## 2. Active / in-flight state

**None.** Every SN / curvilinear / Peierls campaign this agent diagnosed is **merged to
main** (git-verified) — the ERR-026 family (#195/#196/#168/#193/#229/#9), the matvec/sweep
unify (#206), the eager-registry "hang" (#212), the pole-cell NaN (#209), the LD diffusion
limit (ERR-061 #240), the Peierls χ source-index (ERR-063 #257 S10a), the affine flux
typing (#208), the S9 LD boundary-slope verdict.

> Merge-status in memory goes STALE. ALWAYS reconcile a "resume X" / "OPEN #NN" note
> against `git merge-base --is-ancestor <hash> HEAD` and `gh issue view <NN>` before
> acting — never trust a frozen "NOT landed".

**Open, no active work** (pick up only if asked; breadcrumbs in §3): #326 (cylindrical
ordering — adjudicated, remediation open), #123, #128, #132/#100, #129, #170.

## 3. Durable reference (one line each)

- [cylindrical_level_ordering_symmetry_adjudication.md](cylindrical_level_ordering_symmetry_adjudication.md)
  — #326 verdict: the curvilinear MMS is EXACTLY blind to the per-level tie-break; the
  ξ-mirror symmetry adjudicates and says the CLOSURE is broken, not the ordering. Backs L18.
- [coupled_block_operator_numerics.md](coupled_block_operator_numerics.md) — coupled 2×2
  block numerics for `CoupledOperator`: ray/bulk is block-TRIANGULAR (direct), the outer
  `A=M−N` is block-GS; WRAP-first; solve strategy is a FAMILY keyed by spectral character.
- [starting_direction_metric_gauge_derivation.md](starting_direction_metric_gauge_derivation.md)
  — the ψ½ block metric `G_sd` is GAUGE-FREE (any SPD) because `apply_transpose == Aᵀ`;
  `G_sd=0` is the one forbidden value. Dense unit-vector-probe recipe. Backs L17.
- [radial_characteristic_carrier_level_position_key.md](radial_characteristic_carrier_level_position_key.md)
  — the ψ½ direct-march `cells_view(p_idx)` is NOT a p_idx-vs-level bug; carve keys by p_idx.
- [glob_vs_gl_spherical_quadrature_study.md](glob_vs_gl_spherical_quadrature_study.md) —
  Gauss-Lobatto tracks GL at ~1.2× error at resolved N (affordable) but is not a drop-in
  (μ=−1 node ⟹ τ_0=0 ⟹ must be straight-char). Backs L16.
- [curvilinear_inverse_seed_taxonomy.md](curvilinear_inverse_seed_taxonomy.md) — is
  curvilinear `(L+C).solve` an honest SweepOperator? slab yes, cylinder
  quadrature-dependent (product lag is foldable), sphere fixed by route (a). Backs L14.
- [cp_matrix_density_and_sphere_conservation.md](cp_matrix_density_and_sphere_conservation.md)
  — CP [P] is structurally DENSE; bonus OPEN DEFECT: spherical CP breaks row-sum=1 above τ~3.
- [issue_208_flux_displacement_residual_typing_debug_value.md](issue_208_flux_displacement_residual_typing_debug_value.md)
  — the convergence-diagnostic catalogue: FluxDisplacement answers WHERE convergence lags,
  AngularResidual answers WHERE the equation is unsatisfied. Backs L11.
- [curvilinear_tau_clamp_vs_pole_floor.md](curvilinear_tau_clamp_vs_pole_floor.md) — the
  sphere pole-CELL closure is O(h) at r→0, invisible to the volume-weighted L2 gate
  (WONTFIX + an L∞/pole gate). Backs L12.
- [sn_space_angle_discretization_coupling.md](sn_space_angle_discretization_coupling.md) —
  space & angle are separable in Cartesian, genuinely coupled in curvilinear via the M-M thread.
- Open rank-N / Peierls breadcrumbs: [direction_n_quadrature_baseline.md](direction_n_quadrature_baseline.md)
  (#123) · [frame_5_qmc_quadrature.md](frame_5_qmc_quadrature.md) (#128) ·
  [issue_100_class_b_mr_mg.md](issue_100_class_b_mr_mg.md) +
  [issue_132_augmented_nystrom.md](issue_132_augmented_nystrom.md) +
  [issue_132_cylinder_hebert.md](issue_132_cylinder_hebert.md) (#132/#100) ·
  [issue_129_planar_limit.md](issue_129_planar_limit.md) (#129) ·
  [peierls_greens_variant_alpha_decision.md](peierls_greens_variant_alpha_decision.md) +
  [peierls_greens_phase1_closeout.md](peierls_greens_phase1_closeout.md) (Variant-α chain).
