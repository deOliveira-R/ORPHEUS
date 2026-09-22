# Numerics Investigator — Memory Index

One line per entry, a hook and never the content. Behavioural lessons live in `lessons.md`
(read FIRST each dispatch); war stories are cold in `_archive/`; campaign detail lives in the
topic files this index points to.

## 1. Lessons (read first)

- [lessons.md](lessons.md) — **read the file, not this line.** Its top is a SPINE of eight
  meta-lessons (M1 what KIND of question · M2 which LIMIT the claim is exact in · M3 converged
  but wrong is the solver · M4 the degenerate fixture · M5 my instrument lied first · M6 two of
  my measurements contradict · M7 a published claim carries its scope · M8 the reference is the
  first suspect), then L1–L27, whose numbers are STABLE identifiers cited from
  `numerical-bug-signatures` and from a diagnostic script.
- [_archive/](_archive/) — cold: the measured war stories, one file per campaign thread. Never
  loaded; opened when a number needs checking.

## 2. Active / in-flight state

**None.** Every campaign this agent diagnosed is merged; verified against git and `gh`
2026-09-21, including the four whose memory notes were stale (#326, #341, #344 and the
#319/#235 Phase 0 all CLOSED or LANDED, with their gates promoted, not "awaiting promotion").

> Merge status in memory goes STALE. ALWAYS reconcile any "open / owed / awaiting" note against
> `git merge-base --is-ancestor <hash> HEAD` and `gh issue view <NN>` before acting.

**Open, no active work** (pick up only if asked; breadcrumbs in §3): **#343** (the octant sweep
ORDER is an unowned rate lever — the successor #341 left behind), **#319/#235** (later phases of
the cylindrical angular closure), **#200** (block-inverse preconditioner, the seed-lag exit),
#123, #128, #129, #132/#100, #170.

## 3. Durable reference (one line each)

⚠ Five records below live in `scratch/` and are **UNTRACKED** (`git ls-files` 2026-09-21):
they exist on this machine only. The `_archive/` copies of their lessons are the durable half.

- `_archive/issue_344_singularity_kernel_and_gauge.md` — #344: the reflective-box kernel,
  measured then closed-form, and the exact gauge. Backs L23/L24.
- `_archive/curvilinear_tau_ld_and_gram_ownership.md` — τ keeps its arity; the adjoint/Gram
  ownership table. Backs L25/L26.
- `_archive/rate_spectrum_certificate_and_angular_axis.md` — G-S-vs-Jacobi mechanism, the
  refuted outer certificate, the flux-dip axis. Backs L19/L20/L21.
- `_archive/curvilinear_seed_metric_and_ordering.md` — seed taxonomy, GLob study, ψ½ metric,
  level ordering. Backs L14/L16/L17/L18.
- `_archive/krylov_composite_restart_and_stale_references.md` — the grown-composite `restart`
  truncation and the stale-frozen-reference triage. Backs L15/L22.
- `_archive/rcond_threshold_rederivation.md` — the `pinv` cutoff re-derivation. Backs L27.
- [issue_341_boundary_gs_rate.md](issue_341_boundary_gs_rate.md) — DD's undamped `−1` channels
  void Varga; `ndim` refuted; the lever is the octant ORDER (#343).
- [cylindrical_level_ordering_symmetry_adjudication.md](cylindrical_level_ordering_symmetry_adjudication.md)
  — the MMS is exactly blind to the tie-break; ξ-mirror says the CLOSURE is broken.
- [curvilinear_inverse_seed_taxonomy.md](curvilinear_inverse_seed_taxonomy.md) — is curvilinear
  `(L+C).solve` an honest SweepOperator? per (geometry × quadrature).
- [starting_direction_metric_gauge_derivation.md](starting_direction_metric_gauge_derivation.md)
  — the ψ½ block metric is gauge-free (any SPD); `G=0` is the one forbidden value.
- [glob_vs_gl_spherical_quadrature_study.md](glob_vs_gl_spherical_quadrature_study.md) —
  Gauss-Lobatto tracks GL at ~1.2× error but is not a drop-in (`μ=−1` ⟹ straight char).
- [issue_319_flux_dip_discriminator.md](issue_319_flux_dip_discriminator.md) — thickness does
  not split the τ schemes; `h→0` does. β is a sphere-only instrument.
- [curvilinear_tau_clamp_vs_pole_floor.md](curvilinear_tau_clamp_vs_pole_floor.md) — the sphere
  pole-cell closure is O(h), invisible to the volume-weighted L2 gate (WONTFIX + an L∞ gate).
- [issue_208_flux_displacement_residual_typing_debug_value.md](issue_208_flux_displacement_residual_typing_debug_value.md)
  — the convergence-diagnostic catalogue (the typed surfaces have since moved; see L11).
- [coupled_block_operator_numerics.md](coupled_block_operator_numerics.md) — ray/bulk is
  block-TRIANGULAR, the outer `A=M−N` block-G-S; the solve strategy is a FAMILY.
- [radial_characteristic_carrier_level_position_key.md](radial_characteristic_carrier_level_position_key.md)
  — the ψ½ direct march is not a p_idx-vs-level bug; carve keys by p_idx.
- [cp_matrix_density_and_sphere_conservation.md](cp_matrix_density_and_sphere_conservation.md)
  — CP `[P]` is structurally DENSE; OPEN defect: spherical CP breaks row-sum=1 above τ~3.
- [sn_space_angle_discretization_coupling.md](sn_space_angle_discretization_coupling.md) —
  space and angle separate in Cartesian, couple in curvilinear via the M-M thread.
- [cyl_matvec_twin_path_signatures.md](cyl_matvec_twin_path_signatures.md) — the two cylinder
  matvec bugs flat ψ hid for months. Backs L6.
- [krylov_restart_truncation_bug.md](krylov_restart_truncation_bug.md) ·
  [sn_keff_hang_was_eager_registry.md](sn_keff_hang_was_eager_registry.md) ·
  [issue_240_d5b_s3_diffusion_limit.md](issue_240_d5b_s3_diffusion_limit.md) ·
  [atalay_r099_paper_floor_2026_05_03.md](atalay_r099_paper_floor_2026_05_03.md) — the founding
  cases behind L10, L9, L7, L5.
- Open rank-N / Peierls breadcrumbs: [direction_n_quadrature_baseline.md](direction_n_quadrature_baseline.md)
  (#123) · [frame_5_qmc_quadrature.md](frame_5_qmc_quadrature.md) (#128) ·
  [issue_100_class_b_mr_mg.md](issue_100_class_b_mr_mg.md) +
  [issue_132_augmented_nystrom.md](issue_132_augmented_nystrom.md) +
  [issue_132_cylinder_hebert.md](issue_132_cylinder_hebert.md) (#132/#100) ·
  [issue_129_planar_limit.md](issue_129_planar_limit.md) (#129) ·
  [peierls_greens_variant_alpha_decision.md](peierls_greens_variant_alpha_decision.md) +
  [peierls_greens_phase1_closeout.md](peierls_greens_phase1_closeout.md) (Variant-α chain).
