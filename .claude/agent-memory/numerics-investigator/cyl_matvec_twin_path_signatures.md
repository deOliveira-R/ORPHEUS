---
name: cyl-matvec-twin-path-signatures
description: Two reusable cylinder SN matvec bug signatures (legacy code RETIRED; signatures durable) — (1) bool-mask scatter lhs[:,mask] orders by True-position not source-column order = ERR-049 index swap; (2) decoder "analytical extension" at outer-cell inward ordinates diverges O(h) from unified WDD matvec. Both hidden by flat-ψ-only coverage.
metadata:
  type: project
---

# Cylinder SN matvec — two twin-path divergence signatures (Issue #197 / #206)

**Status:** the LEGACY code both signatures lived in is RETIRED on origin/main —
`transport_operator_matvec_cylindrical` and `solution_to_angular_flux_cylindrical`
no longer exist in `orpheus/sn/operator.py` (unified WDD matvec landed).
ERR-049 (Signature 1) is in `tests/l0_error_catalog.md`; the promoted gate is
`tests/sn/test_unified_matvec_cylinder.py` (curvilinear subtree). The remaining
cylinder matvec/sweep-unification work is **Issue #206** (open: "Unify
matvec/sweep WDD recurrence machinery — single source of truth for denom/a/b";
`test_unified_cylinder_matches_hand_reference` is its known-bad). This note keeps
the two reusable SIGNATURES; the per-line code refs are historical.

## WHY BOTH HID — flat-ψ-only cylinder matvec coverage (the meta-lesson)
Every legacy cylinder apply-matvec test used flat ψ (homogeneous reflective +
uniform Q → at convergence ψ flat per ordinate). Under flat ψ BOTH signatures
vanish: Signature-1's permutation becomes the identity (all level-internal values
equal); Signature-2's O(h) decoder mismatch → 0 (cell-average = cell-centre in the
continuous limit). Krylov-on-apply converges to a self-consistent fixed point of
the buggy operator — physically acceptable, still linear+conservative — so nothing
looks wrong. **Cylinder matvec MUST be tested against a per-ordinate hand reference
on NON-FLAT ψ** (the sphere had `test_apply_face_fluxes_match_sweep_recurrence_
spherical`; cylinder lacked the analog). This is vv-principles H2/H3: homogeneous
+ conservation are both degenerate to redistribution/routing bugs.

## SIGNATURE 1 — bool-mask scatter orders by True-position, not source-column order (ERR-049)
`lhs[:, ks] = streaming + redistribution + collision` where
`ks = eq_map.unknowns_at_cell_for_mask(i, mask)` returns ks in ASCENDING global
ordinate order (positions 0..N-1 of `table[i][mask]`), BUT the RHS columns are in
`global_in = level_idx[mask_within_level]` order — sorted by **mu_x within the
level**, NOT globally ascending. When `global_in ≠ sorted(global_in)` the
assignment is a silent PERMUTATION.
- **Fingerprint:** O(magnitude) bit-exact gap on non-flat ψ (rel_diff ~1.1), with
  an EXACT permutation pattern between two mu_x sub-blocks of a level (LS4 p=0
  inward: `global_in=[8,9,10,11,0,1,2,3]` ≠ ascending → ordinates 0..3 (μ=−0.408)
  swap with 8..11 (μ=−0.816); outward bit-exact because LS4 p=0 outward IS
  ascending). Worse for LS6+ (both directions scrambled), LS8+ increasingly so.
- **Sphere immune:** sphere masks on a flat boolean `mu_x < 0` directly on the
  global axis → naturally ascending, no level_idx reordering. The bug is
  cylinder-specific (only cylinder has the level structure that diverges from
  global ascending). H1 "sphere shares the bug" REFUTED.
- **Failure mode #2 (variable/index swap).** Lesson: bool-mask scatter
  `arr[:, mask]` routes by mask-True-position; when source columns come from a
  per-level GATHER, that order may not equal the mask True-position order. ALWAYS
  scatter with explicit fancy-index `arr[:, explicit_indices] = values`.

## SIGNATURE 2 — decoder "analytical extension" diverges O(h) from unified WDD matvec
The cylindrical solution→angular-flux DECODER applied an "analytical extension" at
outer-cell inward ordinates:
`fi[n,:,-1,0] = fi[ref_x[n],:,-1,0]` (reflective-partner OUTWARD cell-centre ψ for
inward n). This fills the cell-CENTRE inward ψ with the outward partner's
cell-centre value — which equals the inflow FACE value only in the continuous
limit. The unified WDD matvec independently computes the outer-cell inward
cell-AVERAGE via its own inward sweep from the BC-traced inflow face. At finite h
the cell-average ≠ outward-partner cell-average by O(h).
- **Fingerprint:** twin-path k_eff disagreement scales O(h):
  `nx=10→2.3e-2, 20→8.5e-3, 40→4.1e-3, 80→2.0e-3`. Matvec residual is machine-zero
  (~6e-13) at the ACTUAL SI ψ but EXPLODES to ~9.7e-2 at the decoder-VIEW ψ
  (replace SI's inward ψ at the outer cell with the reflective-partner outward ψ).
  That contrast PINS the bug to the decoder, not the matvec algebra (the matvec IS
  the sweep⁻¹ on real ψ). The outer-cell rel_diff scales O(h):
  `0.20→0.07→0.03→0.017` across nx∈{10,20,40,80}.
- **The decoder was consistent with the LEGACY symmetric-closure matvec but WRONG
  for the unified WDD matvec** — a twin-path mismatch, not a standalone defect. Fix
  family (preferred): add the outer-cell-inward unknowns to the eq_map so GMRES
  treats them as genuine unknowns and the discrete operator is consistent by
  construction (same shape as the Phase F spherical-pole-IC fix — treat the
  boundary state as a true unknown, NOT a value-based extension). This is the
  Cardinal-Rule-2 / Pattern-2 single-source-of-truth resolution that #206 lands.

## Cross-refs
- ERR-049 is the convention-drift sibling of ERR-050-class operator-algebra↔sweep
  mismatches; see [[r1_step_e_invertible_solve_w_bridge]] for the related
  per-ordinate-vs-iso `/W` convention-drift family.
- Both signatures are instances of the curvilinear-sweep twin-path divergence
  (numerical-bug-signatures Signature 1; ERR-026 manifestation family). The cure is
  Pattern 2: one `SNCellOperator` consumed by BOTH the sweep and the matvec, so the
  paths cannot drift — the architectural target #206 pursues.
