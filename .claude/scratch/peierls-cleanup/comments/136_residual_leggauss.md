# Issue #136 close-out — residual leggauss in slab T_specular + Sanchez K_bc reference (post-#138 docs cleanup)

**Status:** CLOSED — two intentional residual `numpy.polynomial.legendre.leggauss` consumers remain in `orpheus.derivations.peierls_geometry` after the visibility-cone rollout (L3 `e8ad150`): `compute_T_specular_slab` (slab geometric immunity exemption per plan §1B step 3) and `_compute_K_bc_sanchez` (verification-tier reference, out of original plan scope). Both are tracked here as consistency-only contract migrations; no algorithmic change required.
**Date:** 2026-04-30
**Why this comment exists:** Per the post-#138 docs cleanup directive, the §22.9 quadrature-rollout audit in `peierls_nystrom.rst:9910–10216` is being relocated to its originating GitHub issues. This comment carries the **#136-specific rows** (slab T_specular and Sanchez K_bc reference), the geometric-immunity rationale for the slab exemption, and the "out of plan scope" rationale for the Sanchez reference.

<!-- COMMIT-HASH -->

## Per-primitive landing rows for Issue #136

Verbatim from the §22.9 per-primitive landing map (`peierls_nystrom.rst:9939–10026`):

| Primitive | Recipe used | Originating commit | Plan phase |
|---|---|---|---|
| `compute_T_specular_slab` | raw plain GL (geometric immunity — Issue #136) | n/a (plan-exempt §1B step 3) | 1B-skip |
| `_compute_K_bc_sanchez` (verification-tier reference) | raw plain GL (Issue #136) | n/a (out of plan scope) | n/a |

## Residual-consumer rationale (verbatim from §22.9)

Three call sites in `orpheus.derivations.peierls_geometry` remain on raw `numpy.polynomial.legendre.leggauss` after L3 ships; the two pertinent to this issue:

> - `compute_T_specular_slab` (line 2320) — slab exemption per plan §1B step 3 (the integrand $e^{-\tau_{\rm tot}/\mu}$ is geometrically immune at $\mu \to 0$, so plain GL is already spectral). Tracked for consistency-only contract migration in Issue #136.
> - `_compute_K_bc_sanchez` (line 2639) — Sanchez 1986 reference K_bc closed-form, used as a verification benchmark for the cylinder specular closure. Out of original plan scope; tracked in Issue #136 alongside the slab.
>
> The legacy `gl_float` and `gl_nodes_weights` primitives at module-top level remain only as long as those three sites consume them; deletion blocks on Issues #135 and #136 landing.

## Why slab is geometrically immune

The slab specular-transmission integrand is $e^{-\tau_{\rm tot}/\mu}$. The visibility-cone kink that motivated the rollout (chord-impact-parameter discontinuities at shell radii in cylinder/sphere geometry) does not exist in slab geometry — the slab has no shell structure. Consequently plain Gauss-Legendre is already spectral on this integrand, and the rollout plan explicitly skipped slab in §1B step 3:

> *(b)* Phase 4 `compute_T_specular_*` integrals use the new quadrature; N=4 sphere/cyl overshoot magnitude reduced
>
> **DONE** — sphere and cylinder routed through `chord_quadrature` (commit `c4cda20`); slab remains exempt (geometric immunity); the high-$N$ overshoot is now bounded by the structural rank-gating documented at `test_specular_*` (rank-3 within 0.5%, rank-4 within 0.5% gate, `UserWarning` at $N \ge 4$).

Migration to the `Quadrature1D` contract for slab is therefore **consistency-only** — the algorithmic content does not change; the call site would simply consume `gauss_legendre()` (returning a `Quadrature1D`) instead of raw `gl_nodes_weights` tuples.

## Why the Sanchez K_bc reference is out of plan scope

`_compute_K_bc_sanchez` is the closed-form K_bc from Sanchez 1986 used as a verification benchmark for the cylinder specular closure. It sits in the verification tier, not the production tier — the plan was scoped at production primitives and explicitly did not enumerate verification-only references. The same consistency-only argument applies: contract migration would not change the result, just standardise the quadrature interface.

## Production decision

Both call sites remain on raw `numpy.polynomial.legendre.leggauss` — intentional, gated only on consistency cleanup. Deletion of the legacy `gl_float` / `gl_nodes_weights` primitives at the module top of `peierls_geometry.py` blocks on this issue plus Issue #135 (`build_volume_kernel` non-adaptive) both landing.

## Cross-links

- Sphinx (post-cleanup): `docs/theory/peierls_nystrom.rst` §22.9 stub retains the residual-consumer bullet pointing here. Parent recipe abstractions derived in §22.7 (`_section-22-7-visibility-cone:`) — for context on why slab is exempt from the visibility-cone substitution.
- Parent rollout: see Issue #133 close-out for the full §22.9 audit and the acceptance-criterion (b) outcome that pinned the slab exemption.
- Sibling deferrals: Issue #134 (CP-module Phase 1F chord migration); Issue #135 (`build_volume_kernel` non-adaptive).
- Code: `orpheus.derivations.peierls_geometry.compute_T_specular_slab` (line ~2320), `_compute_K_bc_sanchez` (line ~2639), residual primitives `gl_float` / `gl_nodes_weights`. Replacement target: `orpheus.derivations._quadrature.gauss_legendre` constructor.
- Citation: Sanchez 1986 (cylinder specular closure K_bc closed-form reference).
