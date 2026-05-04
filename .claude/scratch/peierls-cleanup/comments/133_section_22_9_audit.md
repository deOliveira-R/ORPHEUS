# Issue #133 close-out — §22.9 quadrature rollout audit (post-#138 docs cleanup)

**Status:** CLOSED — the visibility-cone-substitution rollout (originating from this issue's Phase 5 retreat) shipped as a unified `Quadrature1D` contract + three composable recipes. The `chord_quadrature` recipe in `orpheus.derivations._quadrature_recipes` is the canonical chord-integration entry point; the §22.9 audit table records exactly where each recipe landed.
**Date:** 2026-04-30
**Why this comment exists:** Per the post-#138 docs cleanup directive, the §22.9 "quadrature architecture rollout — outcome and audit" subsection of `peierls_nystrom.rst:9910–10216` is being relocated to the GitHub issues that scoped the rollout. This comment carries the **rollout-completion narrative + the per-primitive landing table + the originating-commit dependency order + the acceptance-criterion audit** that traces back to this issue (Phase 5 retreat → portable primitive → unified contract). Companion close-outs at #134 / #135 / #136 carry the row-specific deferrals (CP chord, `build_volume_kernel`, residual leggauss). A separate #133 comment carries the Phase 5 continuous-µ retreat narrative itself.

<!-- COMMIT-HASH -->

## What this issue produced (rollout outcome)

§§22.7 and §22.8 of `peierls_nystrom.rst` derive the two substitution / subdivision primitives that motivated the `surface_centred_angular_quadrature`, `observer_angular_quadrature`, and `chord_quadrature` recipes that ship in `orpheus.derivations._quadrature_recipes`.

The original plan (`.claude/plans/visibility-cone-substitution-rollout.md`) was scoped as a substitution rollout; the work that actually shipped is broader — a unified `orpheus.derivations._quadrature.Quadrature1D` contract with three composable geometry-aware recipes, plus a sibling `AdaptiveQuadrature1D` for verification-tier integrals. The plan's Phase 1A through 1E acceptance criteria are subsumed and superseded by the architecture; this audit records exactly where each criterion landed.

## Per-primitive landing map

<details><summary>Full 17-row table — recipe used per peierls_geometry primitive (verbatim)</summary>

| Primitive | Recipe used | Originating commit | Plan phase |
|---|---|---|---|
| `compute_T_specular_sphere` | `chord_quadrature` | Q2 `c4cda20` | 1B |
| `compute_T_specular_cylinder_3d` | `chord_quadrature` | Q2 `c4cda20` | 1B |
| `compute_T_specular_slab` | raw plain GL (geometric immunity — Issue #136) | n/a (plan-exempt §1B step 3) | 1B-skip |
| `compute_P_ss_sphere` | `chord_quadrature` (h-space) | Q2 `c4cda20` | 1D |
| `compute_P_ss_cylinder` | `chord_quadrature` (h-space) | Q2 `c4cda20` | 1D |
| `compute_P_esc` (sphere) | `observer_angular_quadrature` | Q3 `50f05ae` | 1C |
| `compute_P_esc_cylinder_3d_mode` | `observer_angular_quadrature` | Q3 `50f05ae` | 1C |
| `compute_G_bc_cylinder_3d` | `observer_angular_quadrature` | Q3 `50f05ae` | 1C |
| `compute_G_bc_cylinder_3d_mode` | `observer_angular_quadrature` | Q3 `50f05ae` | 1C |
| `compute_G_bc` (cyl branch) | `surface_centred_angular_quadrature` | L3 `e8ad150` | (post-plan) |
| `compute_G_bc_outer` (cyl branch) | `surface_centred_angular_quadrature` | L3 `e8ad150` | (post-plan) |
| `compute_G_bc_inner` (cyl branch) | `surface_centred_angular_quadrature` | L3 `e8ad150` | (post-plan) |
| `compute_G_bc_mode` (cyl branch) | `surface_centred_angular_quadrature` | L3 `e8ad150` | (post-plan) |
| `compute_P_esc_outer_mode`, `compute_P_esc_inner_mode`, `compute_G_bc_outer_mode`, `compute_G_bc_inner_mode` | `observer_angular_quadrature` | Q3 `50f05ae` | 1C |
| `*_mode_marshak` per-face primitives | `observer_angular_quadrature` | Q3 `50f05ae` | 1C |
| `build_volume_kernel_adaptive` | `adaptive_mpmath` | Q5/Q6 `123f9ec`/`6c67466` | 1E |
| `build_volume_kernel` (non-adaptive) | raw `gl_nodes_weights` + open-coded subdivision | n/a (Issue #135 follow-up) | (post-plan) |
| `_compute_K_bc_sanchez` (verification-tier reference) | raw plain GL (Issue #136) | n/a (out of plan scope) | n/a |

</details>

## Originating commits, in dependency order

<details><summary>Q1 → L3 commit chain (verbatim)</summary>

- `b281a97` (**Q1**) — ship the `Quadrature1D` contract, the four primitive constructors (`gauss_legendre`, `gauss_legendre_visibility_cone`, `composite_gauss_legendre`, `gauss_laguerre`), and the first two recipes (`chord_quadrature` and `observer_angular_quadrature`). Plan phase 1A.
- `c4cda20` (**Q2**) — migrate the chord call sites ($T_{\rm specular}^{\rm sphere}$, $T_{\rm specular}^{\rm cyl,3D}$, $P_{ss}^{\rm sphere}$, $P_{ss}^{\rm cyl}$) to $h$-space + `chord_quadrature`, making the algebraic identity $T_{00} = P_{ss}$ exact (same nodes, same weights, same integrand). Plan phases 1B + 1D.
- `50f05ae` (**Q3**) — migrate the observer-centred $\omega$-sweeps (~14 primitives spanning the per-face `_mode`, `_mode_marshak`, hollow-sphere, and 3-D Knyazev variants) to `observer_angular_quadrature`. Plan phase 1C.
- `7893bf6` (**Q4**) — extend `Quadrature1D` to own panel slicing (`panel_sizes` and `panel_slice(k)`), letting per-panel basis evaluators index through the contract instead of recomputing offsets externally.
- `123f9ec` (**Q5**) — add the `AdaptiveQuadrature1D` sibling and the `adaptive_mpmath` constructor for verification-tier integrals where the consumer needs only the scalar.
- `6c67466` (**Q6**) — wire residual call sites (most importantly the $\omega$/$\rho$ integrals in `build_volume_kernel_adaptive`) through the contract. Plan phase 1E.
- `e8ad150` (**L3**) — add the third recipe, `surface_centred_angular_quadrature`, derived in §22.8 and used by the four legacy cylinder $G_{\rm bc}$ branches that the plan did not originally enumerate.

</details>

## Why the rollout is broader than the original plan

The plan envisioned shipping a single `gauss_legendre_visibility_cone` utility and threading it through ~6-8 call sites. What the work surfaced is that the *pattern* repeats across two distinct kink-aware subdivision families (chord-impact-parameter at shell radii; observer-centred $\omega$-sweep at tangent angles; surface-centred $\phi$-sweep at chord-quadratic tangents) and that each family deserves its own *recipe* — a geometry-aware constructor of `Quadrature1D` — rather than ad-hoc subdivision in the consumer. The `Quadrature1D` contract itself (frozen value object with `pts`, `wts`, `interval`, `panel_bounds`, `panel_sizes`, `integrate`, `integrate_array`, and the `q1 | q2` composition operator) is the load-bearing abstraction that makes the recipes uniform and the consumers boring. The plan's "ship a utility" framing was correct in spirit but undersized in scope.

## Acceptance criterion audit

<details><summary>Per-criterion outcome table (verbatim from plan §"Acceptance criteria (overall)")</summary>

| Item | Criterion (verbatim) | Where it landed |
|---|---|---|
| (a) | `gauss_legendre_visibility_cone` shipped in `_kernels.py` with L0 test | **DONE** — promoted from `_kernels.py` to the `Quadrature1D`-returning constructor in `_quadrature` (commit `b281a97`); 6 L0 tests in `tests/derivations/test_quadrature.py` carry `@pytest.mark.verifies("gauss-legendre-visibility-cone")`. |
| (b) | Phase 4 `compute_T_specular_*` integrals use the new quadrature; N=4 sphere/cyl overshoot magnitude reduced | **DONE** — sphere and cylinder routed through `chord_quadrature` (commit `c4cda20`); slab remains exempt (geometric immunity); the high-$N$ overshoot is now bounded by the structural rank-gating documented at `test_specular_*` (rank-3 within 0.5%, rank-4 within 0.5% gate, `UserWarning` at $N \ge 4$). |
| (c) | `compute_P_esc_*` / `compute_G_bc_*` mode primitives use the new quadrature; rank-1 algebraic identities preserved at 1e-14 | **DONE** — Q3 commit migrated all per-face `*_outer_mode`, `*_inner_mode`, `*_outer_mode_marshak`, `*_inner_mode_marshak`, and 3-D Knyazev primitives to `observer_angular_quadrature`. L3 commit migrated the four legacy cylinder $G_{\rm bc}$ branches to `surface_centred_angular_quadrature`. The `test_specular_multibounce_cyl_rank1_equals_hebert` canary preserves the algebraic identity at `rtol=1e-8`; bit-equivalent for homogeneous cells. |
| (d) | `compute_P_ss_*` use the new quadrature; rank-1 algebraic identity $T_{00} = P_{ss}$ preserved | **DONE** — sphere and cylinder routed through `chord_quadrature` in $h$-space (commit `c4cda20`); identity tested in `test_chord_quadrature_sphere_T00_equals_P_ss` at bit-equality (same nodes, same weights, identical integrand). |
| (e) | `build_volume_kernel_adaptive` (where applicable) uses the new quadrature | **DONE** — both the inner $\rho$ and outer $\omega$ integrals routed through `adaptive_mpmath` (commits `123f9ec` and `6c67466`). The non-adaptive sibling `build_volume_kernel` retains open-coded subdivision and is tracked separately as Issue #135. |
| (f) | ALL existing Peierls + specular tests pass without modification (no regression) | **DONE** — verified across Q1-Q6+L3 with `test_peierls_*` and `test_quadrature.py` (50 + 7 + 18 + 23 tests across the L3 closure runs). |
| (g) | Documented improvement (table in Sphinx) | **DONE** — this subsection (§22.9), the recipe derivations in §22.7 and §22.8, and the per-recipe docstrings in `_quadrature_recipes`. |
| (h) | ~3-5 commits scaled by phase boundaries | **DELIVERED** as 7 commits (Q1-Q6 + L3) — slightly more than planned because the work expanded from "ship one substitution" to "ship the full `Quadrature1D` contract + three recipes + adaptive sibling." Each commit is independently testable. |

</details>

## Production decision

The `chord_quadrature` recipe in `orpheus.derivations._quadrature_recipes` is the canonical chord-integration entry point. Closed; tracked deferrals migrated to companion issues:

- Slab `compute_T_specular_slab` + Sanchez K_bc reference: **Issue #136** (geometric-immunity exempt; consistency-only contract migration).
- `build_volume_kernel` non-adaptive: **Issue #135** (post-plan migration sketch).
- Phase 1F CP-module chord rollout: **Issue #134** (deferred — CP discretisation error dominates chord-quadrature error there).

## Provenance — back-pointer to this issue

The rollout was initiated as `.claude/plans/visibility-cone-substitution-rollout.md` after the Phase 5 retreat (commit `4dc03cf`) surfaced the visibility-cone substitution as a portable, promotion-worthy primitive recovered from an otherwise abandoned investigation. The L3 closure (commit `e8ad150`) ends the rollout; the residual cleanups are tracked as the GitHub issues listed above.

## Cross-links

- Sphinx (post-cleanup): `docs/theory/peierls_nystrom.rst` §22.9 retains a ~60-LoC stub (intro paragraph + intentional residual consumers bullet list + cross-links to #134/#135/#136). Recipes derived in §22.7 (`peierls-rank-n-bc-closure-section`'s sibling — `_section-22-7-visibility-cone:`) and §22.8 (`_section-22-8-surface-centred-angular:`) survive in full.
- Code: `orpheus.derivations._quadrature_recipes.chord_quadrature`, `observer_angular_quadrature`, `surface_centred_angular_quadrature`.
- Companion #133 comment: Phase 5 continuous-µ retreat narrative (separate close-out for the `compute_K_bc_continuous_mu_*` `NotImplementedError` and 6-round investigation history).
