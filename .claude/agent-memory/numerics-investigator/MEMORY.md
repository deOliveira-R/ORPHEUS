# Numerics Investigator — Memory Index

Slim index. Behavioral lessons live in `lessons.md` (read FIRST each
session). This index holds only (2) live in-flight state and (3) durable
diagnostic references. Bug-campaign play-by-play lives in the topic
`*.md` files this index points to — do NOT inline it here.

## 1. Lessons (read first)

- [lessons.md](lessons.md) — 13 diagnostic-cascade lessons (L1–L13).
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
  `(slice,slice,*idx)`; non-LD tests are blind).

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

- [curvilinear_inverse_seed_taxonomy.md](curvilinear_inverse_seed_taxonomy.md)
  — #226 ruling: is 1-D curvilinear `(L+C).solve` an honest `SweepOperator`? SLAB/CYLINDER
  YES (seedΔ=0, machine residual — cyl α-dome telescopes); SPHERE NO (seed-lagged, residual
  5e5, exact only at FP). Seed-entry map, the seed↔first-ordinate LOCAL CYCLE, feasibility
  (μ=±1 closed → #200 block-inverse precond), Mode-7 MMS-blindness (ansätze ≤linear-in-μ =
  seed's exact regime). Backs L14.
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
