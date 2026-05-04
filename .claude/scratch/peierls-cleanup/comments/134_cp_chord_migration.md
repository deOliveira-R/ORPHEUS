# Issue #134 close-out — Phase 1F CP chord-quadrature migration deferred (post-#138 docs cleanup)

**Status:** CLOSED — Phase 1F (CP-module chord-quadrature migration) was scoped as "optional and lower-priority" by the visibility-cone-substitution rollout plan and remains deferred. CP discretisation error dominates the chord-quadrature error in `orpheus.cp`, so the architectural sharpening yields little observable benefit.
**Date:** 2026-04-30
**Why this comment exists:** Per the post-#138 docs cleanup directive, the §22.9 quadrature-rollout audit in `peierls_nystrom.rst` is being relocated to its originating GitHub issues. This comment carries the **#134-specific deferral rationale** (the "Phase 1F (CP module) — deferred" subsection at `peierls_nystrom.rst:10194–10204`) plus the audit-table row attribution that pinned the deferral.

<!-- COMMIT-HASH -->

## Phase 1F status — verbatim from §22.9

> The plan's Phase 1F was an explicit "optional and lower-priority" rollout of the chord substitution into the `orpheus.cp` module's chord quadrature. CP discretisation error dominates the chord-quadrature error there, so the architectural sharpening yields little observable benefit. Tracked as Issue #134; pick up when a CP regression specifically traces to a chord-quadrature plateau or when a unification sweep across the integral-transport modules is being prepared.

## Originating commit

The **non-CP** chord migration that this issue *would have extended* shipped as commit `c4cda20` (Q2 in the dependency order):

> migrate the chord call sites ($T_{\rm specular}^{\rm sphere}$, $T_{\rm specular}^{\rm cyl,3D}$, $P_{ss}^{\rm sphere}$, $P_{ss}^{\rm cyl}$) to $h$-space + `chord_quadrature`, making the algebraic identity $T_{00} = P_{ss}$ exact (same nodes, same weights, same integrand). Plan phases 1B + 1D.

The audit-table rows pinned to Issue #134 (verbatim from the per-primitive landing map at `peierls_nystrom.rst:9939–10026`):

| Primitive | Recipe used | Originating commit | Plan phase |
|---|---|---|---|
| `compute_P_ss_sphere` | `chord_quadrature` (h-space) | Q2 `c4cda20` | 1D |
| `compute_P_ss_cylinder` | `chord_quadrature` (h-space) | Q2 `c4cda20` | 1D |

These are the **derivations-side** chord migrations that did land. The deferred work is the parallel migration on the **`orpheus.cp` runtime-solver side** (`orpheus.cp.solver` chord quadrature) — phase 1F.

## Acceptance criterion outcome — Phase 1F

The plan's acceptance criterion that this issue tracks (verbatim from §22.9 acceptance audit):

> *(d)* `compute_P_ss_*` use the new quadrature; rank-1 algebraic identity $T_{00} = P_{ss}$ preserved
>
> **DONE** — sphere and cylinder routed through `chord_quadrature` in $h$-space (commit `c4cda20`); identity tested in `test_chord_quadrature_sphere_T00_equals_P_ss` at bit-equality (same nodes, same weights, identical integrand).

That criterion landed for the **derivations-side** chord call sites. Phase 1F's acceptance criterion (which would extend the same migration to the CP-runtime chord quadrature) has no test gate, no commit, and no scheduled work — closed as deferred.

## Production decision

CP discretisation error dominates chord-quadrature error in `orpheus.cp`; sharpening the chord quadrature without first reducing the CP discretisation error is uneconomic. The recipe is available — `orpheus.derivations._quadrature_recipes.chord_quadrature` — and can be wired in if a future CP regression traces specifically to a chord-quadrature plateau, or if a unification sweep across the integral-transport modules is being prepared.

## Cross-links

- Sphinx (post-cleanup): `docs/theory/peierls_nystrom.rst` §22.9 stub retains the deferred-Phase-1F bullet pointing here. Parent recipe derivation in §22.7 (`_section-22-7-visibility-cone:`).
- Parent rollout: see Issue #133 close-out for the full §22.9 audit (rollout outcome, per-primitive landing map, acceptance criterion audit, and the unified `Quadrature1D` contract narrative).
- Sibling deferrals: Issue #135 (`build_volume_kernel` non-adaptive); Issue #136 (residual leggauss in slab `compute_T_specular_slab` and Sanchez K_bc reference).
- Code: `orpheus.derivations._quadrature_recipes.chord_quadrature`; CP-side target is `orpheus.cp.solver` chord quadrature.
