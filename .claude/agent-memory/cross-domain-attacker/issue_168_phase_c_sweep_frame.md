---
name: Issue #168 Phase C sweep-frame matvec
description: Frame-detection precedent — apply matvec restructured as per-ordinate batched sweep over the cell-visit DAG; BC trace law owns the boundary edge; PROMOTES sweep-frame to A.7 sharpening
type: project
---

Issue #168 Phase C structural analysis (2026-05-12). The orchestrator
asked whether the apply matvec's `BoundaryFaceFlux.DDExtrapolation`
(Phase A) and the BC-fed-cell-centres in `solution_to_angular_flux_spherical`
constitute an architectural inconsistency vs. the Wave 0-12 trace-space
+ boundary-law machinery.

**Confirmed inconsistency** (operator.py:533 `bc.apply` consumes
cell-centres at i=N-1; operator.py:701 `DDExtrapolation` ignores BC
entirely — both bypass the §16A.3 affine BC form which requires the
boundary TRACE).

**Best frame: Sweep / wavefront frame composited with functional-
analytic frame.**

- Trigger: `creates_sweep_cycle` ClassVar already exists in
  `BoundaryTraceLaw` and signals exactly the structural property — BCs
  introduce DAG cycles. For matvec (input ψ fixed), the cycle linearises
  to a partner-ordering dependency, not a fixed-point iteration.
- Trigger: `InflowTraceSpace`/`OutflowTraceSpace` exist as first-class
  typed function spaces with per-face directional masks. Ghost-cells
  map to typed trace-space vectors, NOT extrapolated virtual cells.

**Concrete reformulation**:
1. Matvec becomes per-ordinate batched sweep over
   `SNMesh.iter_cell_visits(ordinate_idx=n)`.
2. WDD diamond `ψ_face_out = 2·ψ_cell − ψ_face_in` evaluated per cell.
3. BC trace law applied ONCE at the boundary edge as
   `bc_outer.apply(outflow_trace)` — the outflow trace IS the
   WDD-propagated outflow face, not the input cell-centre.
4. `BoundaryFaceFlux` Protocol + `DDExtrapolation` + `CellCenter`
   retire entirely (Phase A patch subsumed by sweep frame).
5. `PoleAngularClosure` Protocol stays (pole is intrinsic geometry, not BC).

**Elegance payoff** (hits 4 criteria):
- Structure-exposing: the `creates_sweep_cycle` machinery the sweep
  planner already consumes (§15A.2) becomes the matvec's machinery.
- Expressive: matvec becomes "one sweep iteration" semantically.
- Structurally-simpler: 3 defects close by construction in 1 rewrite.
- Algorithmic-advantage: one closure rule covers boundary + interior;
  per-closure (LD, EC, Step) extension becomes a pluggable `E` swap.

**Why this matters for the trigger table**:
The sweep-frame for matvec is a high-value precedent — anytime a
project ships TWO operator paths (sweep vs apply, or solve vs matvec,
or in MoC: track sweep vs matrix-free apply), the sweep frame is
likely the unifying frame if the alternate path implements the same
discrete operator via different storage conventions. Add to reference.md
Part A.7 (Integral equations and special structure):

  TRIGGER: "Two paths to the same discrete operator over different
  storage conventions, with order-degradation at boundaries on one
  path."
  FRAME: Sweep / wavefront — recover the boundary as a DAG edge
  consumed via the trace law, not as a cell-centre algebraic closure.

Add to Part C elegance smells:

  SMELL 16: Structurally distinct paths to the same operator. When
  `apply` and `solve` are supposed to be `L` and `L⁻¹` for the SAME
  operator L, but they operate on different discrete representations
  (cell-centres for one, faces for the other), the operator-correctness
  claim is at risk. The fix is to make both consume the same primary
  representation — typically the faces (cells are derivatives of faces
  under DD).

**Promotion candidate**: Phase C ships → promote to
scripts/validated_sn_sweep_frame_matvec.md with quantitative
order-of-convergence evidence (Phase C Gate 3.1 / 3.2 MMS tests).

Cross-link to:
- elegance_smell_rank_non_monotone.md (Smell 15) — Smell 16 is a sibling.
- variant_alpha_2surface_bie_frame.md — the rank-1/rank-2 boundary-
  scattering operator framing has the same shape (BC trace as the
  off-diagonal block of a resolvent).
