---
name: discretization-scheme-naming-signal
description: NAMING-SIGNAL frame attack on orpheus/sn/spatial/scheme.py — per-name keep/rename verdicts for the model-agnostic advection–reaction spatial-closure interface (consumed soon by a diffusion solver). New smell candidate: frame-leak parameter naming.
metadata:
  type: project
---

# DiscretizationScheme naming-signal assessment (branch feature/sn-space-angle-tier2)

VOCABULARY frame-detection (math settled). The SN within-group
`μ ∂ψ/∂x + Σ_t ψ = Q` IS a linear ADVECTION–REACTION equation; the spatial
closures are generic advection-discretization schemes (Step↔first-order
upwind/donor-cell; DD↔central/Keller box; LD↔DG-P1-upwind). The per-cell
triple `(a, inverse_denom, w)` is generic in (wave-speed `|μ|→a`,
reaction-rate `Σ_t→r`). **A diffusion solver will soon consume this layer
(standalone + consistent-DSA preconditioner)** → interface should read
model-agnostic to a CFD/numerical-PDE reader, to the extent it doesn't lose
precision.

**Why:** the user is framing a docstring + deciding a DEFERRED rename; this
is the durable verdict set so a fresh session doesn't re-derive it.
**How to apply:** when the rename lands, use this table; when the diffusion
consumer is wired, run the name-#1 first-test (does it pass Σ_r/Σ_a?).

## Per-name verdicts (priority order)

1. **`total_xs`/`sigt_cells`/`sig_t` → RENAME `reaction_xs`/`removal_xs`,
   collapse 3 spellings to 1** (HIGHEST). Generic role = zeroth-order
   reaction coefficient `r` (LeVeque 2002 §E; Morton & Mayers). Diffusion
   passes Σ_r/Σ_a here → `total_xs` becomes a LIE under the 2nd consumer.
   Smell #16 shape-2 (one quantity, 3 spellings). FIRST TEST: grep the
   diffusion-DSA call site — if it passes a removal/absorption σ, rename is
   forced. (Which generic spelling = one-grep decision; not yet resolved.)
2. **`cell_average_weight` (long handle for `w`) → RENAME `face_blend_weight`**
   (keep `w` as math symbol). It's the CFD κ-scheme/blended-convection axis
   (van Leer 1979; Hirsch Vol.2 §8): upwind(1)↔central(½). Current name
   names the OUTPUT not the AXIS. **Do NOT use `kappa`** — different
   normalization (3-pt slope vs 2-pt convex weight) = precision-LOSING import.
3. **`streaming`/`s_axes` → KEEP — FRAME-CONFLICT flag.** `|μ|/Δ` is the SN
   streaming coefficient (Lewis & Miller). Diffusion has NO advective
   μ-coefficient (its flux is −D∇φ) → won't populate this slot. Genericizing
   to `advective_coeff` MISDESCRIBES for diffusion AND de-signals for SN.
   The ONE name where SN/CFD genuinely diverge; keep SN-native. Future
   model-agnostic pass must NOT genericize this by reflex.
4. **`a`/`a_attenuation` → KEEP.** Transmission/attenuation multiplier;
   dual-readable (transport + FV/DG local transfer).
5. **`inverse_denom` → KEEP** (`inverse_cell_balance_diagonal` marginally
   sharper, verbose; low priority).
6. **`source_emission`/`cell_average`/`outgoing_face_from_average` → KEEP.**
   Cleanest names in the set; `outgoing_face_from_average` = explicit inverse
   face-reconstruction, direction-clear.
7. **`DiscretizationScheme`/`Base` → KEEP, DEFER.** `SpatialClosure`/
   `CellClosure` is higher-signal ("closure" = precise FV/SN term, LM-1989
   §5.3; under-discriminates vs angular/energy discretizations) but degree-43
   hub → churn not worth modest gain. Put "spatial closure (Keller box /
   DG-upwind / donor-cell)" in the DOCSTRING instead.
8. **`affine_scan_coefficients`, `cell_kernel_batch`/`residual_kernel_batch`,
   `is_affine_scannable`/`is_linear`/`is_positivity_preserving`,
   `Q`/`source`/`psi_*`/`face_*` → KEEP.** Nit: `psi_bar`/`cell_avg`/
   `cell_average_flux` = 3 spellings of the cell average (Smell-16 shape-2
   lexical drift); pick one per layer + document.

## NEW SMELL CANDIDATE — "frame-leak parameter naming" (promotion candidate, Part C)

A model-AGNOSTIC interface that names a parameter after ONE consumer's
physics (`total_xs` = "cross section" on an advection–reaction layer that a
diffusion solver will also consume) becomes a LIE the moment the second
consumer passes a different realization through the same slot. The name
asserts a model the interface no longer owns. TELL: a parameter docstring
that says "generic in X" while the parameter is named after a specific X₁.
FIX: name the ROLE (`reaction_xs`), not the realization (`total_xs`).
Distinct from Smell #16 (that's two representations of one quantity; this is
one slot whose NAME over-commits to one of N consumers). Cross-link:
coding-elegance "code should read like the math/domain" — here the domain is
plural (transport ∪ CFD), so the name must read in the INTERSECTION.

## Checked-and-confirmed (not speculation)
- advection–reaction frame: CONFIRMED right (literature this session; the
  file ALREADY documents it at scheme.py:568-570, 827). Highest-signal
  frame for this vocabulary.
- κ-scheme frame: CONFIRMED relevant for `w`, but `kappa` itself is a
  precision-losing name (normalization mismatch) — frame near-miss.
- streaming↔CFD-advective: CONFIRMED divergent (frame conflict), keep SN.

## Files
- /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/spatial/scheme.py
- /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/spatial/diamond.py
- /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/spatial/linear_discontinuous.py
