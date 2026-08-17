# Numerics Investigator — Memory Index

One line per entry. Behavioral lessons live in `lessons.md` (read FIRST each session);
campaign detail lives in the topic files this index points to — never inline it here.

## 1. Lessons (read first)

**⚠ READ `lessons.md` ITSELF — it is the authority; the keywords below are a lookup
table, not a summary.** The spine: never guess, isolate in cascade order, and a single
mesh / flat flux / homogeneous case / two-probe agreement proves nothing.

- [lessons.md](lessons.md) — L1–L24. Keyword index (one line each, `##` headings there):
  - L1 cascade order · L2 curvilinear redistribution is prime suspect · L3 two-quadrature
    signed-error gate · L4 convergence-rate fingerprints · L5 paper-floor vs code-bug
  - L6 curvilinear matvec needs a NON-FLAT reference · L7 per-ordinate moment → GLOBAL
    frame before the angular reduction · L8 our own theory page can be the contamination
  - L9 a "hang" may be fixture cost · L10 diverges-with-refinement + a discarded info-flag
    = unconverged inner (**b** recorded-but-unread flag; **c** all-reflective absorber is
    SI-hard, `Σ_t·n_inner` invariant, G-S-vs-Jacobi INVERTS with dimension)
  - L11 residual `r=Aψ−q`, not the increment · L12 offline error is THE floor only after
    swap + silent control + AMPLIFY · L13 greedy `(Ellipsis,*idx)` under a spectator axis
  - L14 curvilinear seed-lag is QUADRATURE-dependent · L15 sweep ANGULAR N (not h) for a
    closure re-pose; resize Krylov `restart` from the composite · L16 cross-quadrature
    comparison needs a scheme-faithful driver + contamination guard
  - L17 a STATE DOF's metric ≠ its angular weight; `T=Aᵀ` ⟹ gauge-free, `G=0` forbidden
  - L18 adjudicate a LABELING degeneracy with the operator's own SYMMETRY GROUP
  - L19 a RATE question is a SPECTRUM question (eigen-solve `M⁻¹N`); check the comparison
    THEOREM before hunting; enumerate a finite knob; never branch on a correlated variable
  - L20 a RESIDUAL cannot gate an EIGENVALUE contract — measure the TRANSFER GAIN first
  - L21 an ANGULAR claim separates by `h→0`, not by the parameter it is named after
  - L22 a STALE frozen reference: relative error before nulp; a re-baseline's radius is
    FROZEN REFERENCES, not `.npz` files
  - L23 a SINGULARITY is a two-question object (refuse the either/or); dense SVD through
    the production builders settles what ARPACK bounds; the METRIC picks the remedy;
    fit → predict → test OFF-SAMPLE → SWAP THE SCHEME. **+addendum (disposition):** change
    the SPLITTING not the mesh; a ladder can be PARITY-SPLIT; enumerate the moment ladder.
    **+addendum II (coherence):** remove the kernel ≥3 ways; `‖MM⁻¹−I‖` full-space is the
    wrong instrument for a forward-substitution (SUBSPACE) inverse
  - L24 a KERNEL is usually a CLOSED-FORM problem — a counting law INDEPENDENT of a
    parameter the operator contains ⟹ the governing equation is combinatorial; sign
    CHARACTERS diagonalise a specular-BC system; a SUM over axes = modes on PLANES; where
    SVD is unaffordable the span check is a PRODUCTION-generated kernel vector + a
    round-off NEGATIVE control; ABSENT vs INVISIBLE share a parity fingerprint; a
    blindness list measured on ONE quadrature is a sample; never ship a DENSE basis

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
ordering — adjudicated, remediation open), **#341** (boundary-G-S rate — mechanism found,
`ndim` refuted, default unchanged; docstring fixes + the octant-order issue owed),
**#344** (ANSWERED 2026-08-14, THREE halves — structure, disposition
(DETERMINISM / solver-selection, exactly removable) and COHERENCE (boundary-G-S IS
a splitting of A: `M−N≡A` bit-exact, kernel-free controls agree at 1e-13, ERR-056
positive control at 1.0); remediation is the user's call, 24 gates awaiting
promotion; see §3),
**#319/#235** (flux-dip experiment RUN and answered — see §3; remediation is the user's
call, 3 gates awaiting promotion), #123, #128, #132/#100, #129, #170.

## 3. Durable reference (one line each)

- [issue_341_boundary_gs_rate.md](issue_341_boundary_gs_rate.md) — #341 verdict: DD's `d−1`
  undamped `−1` channels void Varga ⟹ G-S may lose; `ndim` REFUTED as the predicate; the
  real lever is the unowned octant ORDER (`max_a L_a > Σ_{b≠a} L_b`, 25/25). Backs L19.
- ⭐ `scratch/issue_344_null_space_structure.md` (in-repo, not agent memory) — #344, BOTH
  halves. **Part I (structure):** `ker A = ` tangential slots (`G ≡ 0`, `product`/`lebedev`
  only) `⊕` a real DD underdetermination (`ng·N/4` at d=2, `ng·(N/8)·(2Σn−1)` at d=3,
  `level_symmetric` has ONLY this); LD non-singular on the same box; every `d ≥ 2`
  Cartesian DD reflective k-solve is singular by default. **Part II (disposition):
  DETERMINISM, not a discretization bug** — the deviation is boundary-G-S's oblique
  projector (Jacobi returns the exact trace, 5/5), converges at exactly `O(h)`, and the
  G-orthogonal gauge projection is an EXACT fix (`8.97e-02 → 5.8e-13`).
  **Part III (coherence):** the schedule IS a splitting — `M−N ≡ A` bit-exact 20/20,
  four kernel-free controls agree at `1e-13`, ERR-056 mutation at `1.0`; `M_GS⁻¹` is
  a SUBSPACE inverse. ⛔ the `n_x`-parity split is SOURCE-dependent (aniso: `1.3e-02`
  at even `n_x`). Backs L23 + both addenda.
- ⭐ `scratch/issue_344_kernel_basis.md` (in-repo) — #344's **CLOSED-FORM basis** for
  `ker A`: bulk-zero face sawtooth `(−1)^k φ_a(i_⊥)` ⟹ `Σ_a s_a Y_a(s_{≠a};i_{≠a}) = 0`,
  basis = (group)×(ordinate orbit)×(axis PAIR)×(sign character)×(complementary cell index).
  Both counting laws become theorems; verified vs dense SVD on 13 configs (res `≤1.6e-15`,
  subspace gap `≤2.2e-14`). Gauge projection SHIPPABLE: `11.26 % → 5.8e-13`, bulk
  bit-unchanged, setup `0.04 s` / apply `0.094 ms`, blocked storage (17.6 GiB → 154 MiB).
  ⛔ `p=0` face moment is NOT blind when tangential slots exist. Backs L24.
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
- [issue_319_flux_dip_discriminator.md](issue_319_flux_dip_discriminator.md) — #319 verdict:
  the thickness DECAY RATE does not split M-M τ from plain diamond (0.000 for both); `h→0`
  does (3.2×→204×, λ_opt→1). Reconciles #235 as wrong-regime (`c=0.5`) + high angular order
  (14×@S2 → 0.9×@S8). β is a sphere-only instrument. Backs L21.
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
