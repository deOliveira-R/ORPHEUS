---
name: issue-240-d5b-s2-ubld-wiring
description: D5b-S2 d≥2 UBLD LD kernel wiring — PASS-WITH-NITS; the moment-axis-tail twin + slot-0 convention dup; orphaned docstring; honest scalar-Q kernel (crosswalk stale)
metadata:
  type: project
---

# #240 D5b-S2 — the d≥2 UBLD Linear-Discontinuous kernel wiring

**Verdict: PASS-WITH-NITS** (branch `feature/sn-space-angle-tier2`, reviewed pre-commit).
Closes LD's d≥2 raise so LD runs on the DAG wavefront in 2-D via the verified
`_ubld.py` per-cell primitive (`assemble_ubld`/`per_cell_solve`/`assemble_inflow_axis`).
The reduction-lattice (`spatial_basis_per_axis` ClassVar, `per_axis**d` cell /
`per_axis**(d-1)` face) is a GENUINE reduction, NOT a forked twin — DD/Step stay
byte-identical by construction (default trait =1 → `tail=()`, no length-1 axis;
verified 38 spatial+primitive+MMS-foundation green, protocol+dispatch 38 green,
cartesian_2d collected-0 [file moved]).

## What is ELEGANT here (reinforce — the implementer nailed these)
- **Q_cells handling is HONEST, not shape-sniffing.** The crosswalk PROPOSED "kernel
  accepts BOTH rank-r and rank-(r+1) Q_cells" (a boolean-flag-in-disguise). The CODE
  rejected that: `_ubld_system` UNCONDITIONALLY lifts the scalar avg to slot 0
  (`S_moments[...,0] = Q + zeros(batch)`), full moment-vector deferred to S3. Code is
  MORE elegant than its design record. The crosswalk's "accepts BOTH" is now stale-on-paper.
- **d=1 fast-path cleanly separated.** `cell_kernel_batch`/`residual_kernel_batch` both
  branch `if len(s_axes)==1:` → `_kernel_terms`→`d1_closed_form` (NO `per_cell_solve`);
  d≥2 → dense. d=1 never routes through the dense solve (L16 preserved). `_kernel_terms`
  raise changed NotImplementedError→ValueError, now unreachable-defensive (dispatch is
  the real gate); only internal callers, all len-1-guaranteed.
- **`_CellResidual` d≥2 raises NotImplementedError (S3 deferral) — illegal-state-LOUD.**
  The multi-D matvec needs the 2^d spatial iterate (S3); rather than feed a scalar probe
  that broadcasts wrong, it fails loudly. The kernel-level matvec twin
  (`residual_kernel_batch` at the full solved state) IS verified directly (D5b.1/D5b.4).
- **Closure-consistency (out-of-cell==in-of-next-cell) is the load-bearing convention**
  and is pinned: `_ubld_outgoing_faces` sums o_a∈{0,1} blocks in the SAME transverse-
  Kronecker order `assemble_inflow_axis` consumes; round-trip D5b.1 + matvec-twin D5b.4
  pin it with NON-flat per-axis psi_in + both-direction face reconstruction.

## NITS to fix before commit (actioned by main agent)
1. **MUST-FIX (Cardinal Rule 3 doc regression): orphaned `is_affine_scannable` docstring.**
   The new `spatial_basis_per_axis = 2` + its docstring were inserted BETWEEN
   `is_affine_scannable = True` (linear_discontinuous.py:284) and its original docstring.
   AST-verified: `is_affine_scannable` now has NO following string literal; its docstring
   (the load-bearing "LD rides the DAG wavefront in d≥2, not the scan-march" rationale —
   the EXACT trait this carve relies on) sits as a dead 2nd consecutive literal under
   `spatial_basis_per_axis`. Fix: move the `is_affine_scannable` docstring back under its
   attribute (above the `spatial_basis_per_axis` line).
2. **CONCERN (Pattern 7 + Pattern 3): slot-0 "average moment" convention re-spelled at 6
   sites** across 3 files (`_ubld_system:533`, sweep_graph `_CellSolve:879`,
   loss_rep `_inflow_to_moments:313`/`octant_pass:962`/`_octant_face_cochain:1098`/
   `_edge_outflow:1194`). The Kronecker ordering `[bar,ŷ,x̂,x̂y]` puts the average at index
   0 (pinned in `_ubld.py` + test_ld_ubld_primitive `psi_exact`), but nothing NAMES it.
   Two sites (`capture[...,0]` in both reps) are verbatim twins. Fix: single-source a named
   `AVERAGE_MOMENT = 0` in `_ubld.py`, reference everywhere. Bug habitat: a basis-ordering
   change moves the average off slot 0 → 6 silent miscatches.
3. **CONCERN (twin-storage idiom): `tail = () if n_face_moments==1 else (n_face_moments,)`
   duplicated** in `_MovingFrontier.__init__` (sweep_graph:267) AND
   `FullFieldWavefront._octant_face_cochain` (loss_rep:1084) — the two storage policies
   (storage-B rolling window / storage-A full cochain) for the SAME cochain. Value
   single-sourced (`n_face_moments` threaded from `_n_face_moments`); the value→shape-suffix
   MAPPING is independently spelled in two MODULES. Classic divergence habitat (matches the
   institutional "two storage policies, independently hand-edited" tell). Fix: extract
   `_face_tail(n)->tuple[int,...]` shared, or one buffer-builder owning the notion.
4. **CONCERN (two-spellings-of-one-partition): the multi-moment gate has 2 spellings.**
   `_CellSolve`/`_CellResidual` gate `len(s_axes)>1 and per_axis>1`; loss_rep gates
   `_n_face_moments != 1` (= `per_axis**(ndim-1) != 1` = `per_axis>1 and ndim>1`).
   Coextensive ONLY under the invariant `len(s_axes)==ndim` (holds: walk builds s_axes via
   `zip(str_axes_octant, cell_idx)`, one per spatial axis). STRENGTH today (shared concept);
   smell is the 2nd spelling. Bug habitat: a future partial-s_axes caller (single-axis
   scan-march on a multi-D mesh) breaks coextensivity. LD is the scan-march EXCL for d≥2
   so it can't arise now — CONCERN not VIOLATION.

## Carry-forward for S3/S4 reviews
- S3 wires the 2^d spatial-moment iterate φ̂ → removes the `_CellResidual` d≥2 raise +
  threads the slope SOURCE (Q̂≠0). The slot-0 convention (nit 2) WILL spread further then —
  single-source it NOW.
- S4 widens BoundaryFlux for non-vacuum domain inflow → `_inflow_to_moments`/`_edge_outflow`
  stop being vacuum-zero identities. The `[...,0]` capture extractions become live.
- d=2 numpy↔symbolic A==A entry-wise pin LANDED (`test_d2_assembled_matrices_match_symbolic`,
  catches ERR-060) — the S1-hindsight carry-forward is closed.
