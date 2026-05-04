# Issue #135 close-out — `build_volume_kernel` non-adaptive migration deferred (post-#138 docs cleanup)

**Status:** CLOSED — `build_volume_kernel_adaptive` was migrated to the unified `Quadrature1D` contract via `adaptive_mpmath` (commits Q5 `123f9ec` and Q6 `6c67466`). The non-adaptive sibling `build_volume_kernel` predates the recipe abstraction and was not enumerated in the original visibility-cone rollout plan; it remains on raw `gl_nodes_weights` + open-coded subdivision, intentionally deferred as a residual consumer.
**Date:** 2026-04-30
**Why this comment exists:** Per the post-#138 docs cleanup directive, the §22.9 quadrature-rollout audit in `peierls_nystrom.rst:9910–10216` is being relocated to its originating GitHub issues. This comment carries the **#135-specific row** (the `build_volume_kernel` non-adaptive entry), the originating Q5/Q6 commits (which migrated the *adaptive* sibling), and the acceptance-criterion audit row (e) that scoped the deferral.

<!-- COMMIT-HASH -->

## Per-primitive landing rows for Issue #135

Verbatim from the §22.9 per-primitive landing map (`peierls_nystrom.rst:9939–10026`):

| Primitive | Recipe used | Originating commit | Plan phase |
|---|---|---|---|
| `build_volume_kernel_adaptive` | `adaptive_mpmath` | Q5/Q6 `123f9ec`/`6c67466` | 1E |
| `build_volume_kernel` (non-adaptive) | raw `gl_nodes_weights` + open-coded subdivision | n/a (Issue #135 follow-up) | (post-plan) |

The adaptive sibling was migrated; the non-adaptive variant **is the deferred work this issue tracks**.

## Originating commits (the migration that DID land — adaptive sibling)

Verbatim from the §22.9 originating-commits dependency order:

- `123f9ec` (**Q5**) — add the `AdaptiveQuadrature1D` sibling and the `adaptive_mpmath` constructor for verification-tier integrals where the consumer needs only the scalar.
- `6c67466` (**Q6**) — wire residual call sites (most importantly the $\omega$/$\rho$ integrals in `build_volume_kernel_adaptive`) through the contract. Plan phase 1E.

## Acceptance criterion — Phase 1E outcome

Verbatim from the §22.9 acceptance-criterion audit:

> *(e)* `build_volume_kernel_adaptive` (where applicable) uses the new quadrature
>
> **DONE** — both the inner $\rho$ and outer $\omega$ integrals routed through `adaptive_mpmath` (commits `123f9ec` and `6c67466`). The non-adaptive sibling `build_volume_kernel` retains open-coded subdivision and is tracked separately as Issue #135.

## Residual-consumer status (verbatim from §22.9)

Three call sites in `orpheus.derivations.peierls_geometry` remain on raw `numpy.polynomial.legendre.leggauss` after L3 ships; the row pertinent to this issue:

> `build_volume_kernel` (line 1194, the non-adaptive variant) — predates the recipe abstraction and wasn't enumerated in the plan. Tracked as Issue #135 with a migration sketch.
>
> The legacy `gl_float` and `gl_nodes_weights` primitives at module-top level remain only as long as those three sites consume them; deletion blocks on Issues #135 and #136 landing.

## Migration sketch (the recipe that would replace the open-coded subdivision)

The non-adaptive `build_volume_kernel` consumes raw `gl_nodes_weights` + an open-coded subdivision of the integration domain. The natural replacement is `composite_gauss_legendre` (built on the `Quadrature1D` contract via the `q1 | q2` composition operator) plus, where the integrand has visibility-cone kink structure, `observer_angular_quadrature`. Picking up this work would:

1. Identify each subdivision boundary in the open-coded `build_volume_kernel` body.
2. Replace the manual offset arithmetic with `Quadrature1D.panel_slice(k)` (commit `7893bf6` extended the contract to own panel slicing for exactly this consumer pattern).
3. Verify bit-equality (or convergence-equivalence) against the `_adaptive` reference.
4. Remove the open-coded subdivision; `gl_float` / `gl_nodes_weights` deletion unblocks once #135 + #136 both land.

## Production decision

`build_volume_kernel` non-adaptive remains on raw GL + open-coded subdivision — this is intentional: it is a residual consumer of the legacy `gl_nodes_weights` primitive, kept until the migration sketch above is executed. The adaptive sibling carries the production load for verification-tier integrals.

## Cross-links

- Sphinx (post-cleanup): `docs/theory/peierls_nystrom.rst` §22.9 stub retains the residual-consumer bullet pointing here. Parent recipe abstractions derived in §22.8 (`_section-22-8-surface-centred-angular:`).
- Parent rollout: see Issue #133 close-out for the full §22.9 audit and the unified `Quadrature1D` contract narrative.
- Sibling deferrals: Issue #134 (CP-module Phase 1F chord migration); Issue #136 (residual leggauss in slab `compute_T_specular_slab` and Sanchez K_bc reference).
- Code: `orpheus.derivations.peierls_geometry.build_volume_kernel` (deferred, line ~1194), `build_volume_kernel_adaptive` (migrated). Recipe library: `orpheus.derivations._quadrature_recipes`. Contract: `orpheus.derivations._quadrature.Quadrature1D` / `AdaptiveQuadrature1D`.
