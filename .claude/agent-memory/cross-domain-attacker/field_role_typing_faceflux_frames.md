---
name: field-role-typing-faceflux-frames
description: Durable frame for the SN face-flux storage locus — interior FaceFlux is a primal 1-cochain (DEC); the boundary trace is its restriction ι* to the boundary 1-chain, NOT a distinct quantity; "absorption = identity" IS ι_*∘ι*=id. Verdict: per-face object REJECTED (vectorization); field-on-FaceLayout+views is native (biproduct C¹=C¹_int⊕C¹_∂). Second frame: sparse triangular factorization makes G-S/Jacobi a back-edge split of the reflective coupling on the (octant×face) graph.
metadata:
  type: project
---

# Field-role-typing — FaceFlux storage-locus frame attack

The question: is a face-located, ephemeral, single-role `FaceFlux` (with the boundary
trace as its domain-edge subset) the native decomposition, or does it fight the
problem's frame? Two native frames fire. (Branch namesake; the cochain reading landed.)

## STRUCTURAL FEATURES (the triggers)
1. The quantity lives on **faces (codim-1)**; the cell-average ψ̄ (codim-0) is derived
   from the faces by the diamond closure `ψ̄=(ψ_in+ψ_out)/2`.
2. The face flux is **directional per ordinate**: `sign(Ω·n_face)` selects in/out face
   AND orients the cell DAG AND partitions the inflow/outflow trace — ONE sign function,
   three consumers.
3. The sweep is **forward substitution over an upwind DAG** — inversion of a sparse
   lower-triangular `L_oct` in the upwind cell ordering.
4. The boundary is the **codim-1 locus where the face flux meets the BC**; reflective
   coupling is a permutation `Ω↔Ω'` between boundary faces.
5. The diamond closure is a **local two-term incidence relation** (`ψ_out=2ψ̄−ψ_in`).
6. TWO persistence strata: boundary slots persist across sweep calls; interior slots
   are ephemeral (rebuilt per call).
7. The G-S/Jacobi question is a **splitting of the reflective back-edge set** of the
   face-coupling graph.
8. **Vectorization is load-bearing**: the unit of operation is the whole
   `(N_oct, ng, n_diag)` wavefront batch; a per-cell/per-face Python fold cost a
   10–20× regression.

Ground-truth fact that clinches Frame 1: the boundary↔interior exchange is a LITERAL
identity copy at the domain-edge slot (the sweep seeds `psi_x[…,0,…]=boundary.face_view`
and absorbs `boundary.face_view[:]=psi_x[…,0,…]`). The domain-edge interior face IS the
boundary trace slot. "Absorption = identity" is exactly literal.

## Frame 1 — Discrete exterior calculus / cochains (THE native frame)
Trigger: features 1 + 5 + 4. The angular flux crossing faces, per ordinate Ω, is a
**primal 1-cochain** `ψ⁽¹⁾_Ω ∈ C¹(mesh)`. The cell-average ψ̄ is a 0-cochain; the
diamond closure is the averaging (Hodge-star / projection) map `⋆:C¹→C⁰`. The boundary
trace is the **pullback** `ι*ψ⁽¹⁾` under the inclusion `ι:∂Ω↪Ω` — exactly the
"select the domain-edge faces" operation the sweep does by hand. The single unified
object: ONE face-cochain field on a FaceLayout enumerating ALL faces, with a
distinguished `boundary_view()` = `ι*`. "Absorption = identity" because `ι_*∘ι*=id`
on the boundary 1-chain. **The proposal IS the discrete trace operator.**

Elegance payoff: structure-exposing (manual seed/absorb marshalling becomes `ι*`/`ι_*`;
the two persistence strata become the chain decomposition `C¹=C¹_int⊕C¹_∂`, not an
ad-hoc storage policy), expressive (ONE field type replaces typed-boundary + raw-numpy),
structurally-simpler (collapses two storage representations into one). Algorithmic-
advantage does NOT fire — the wavefront DAG already gives the optimal schedule; this is
a representation win, not a speed win. State that plainly; do not invent a perf payoff.

First test (must fail-able): build a full-face FaceLayout + FaceFlux; assert that seeding
via `boundary_view().copy_from(inflow)` and reading it back after the sweep reproduces
the current boundary buffer BIT-FOR-BIT. Fails if there is any reindexing/transpose
between the interior-edge face ordering and the trace layout ordering (a real risk).

## Frame 2 — Sparse triangular factorization / transfer-matrix (the G-S/Jacobi illuminator)
Trigger: feature 3 + feature 7. Per octant, the streaming operator is a sparse
**lower-triangular** `L_oct` in upwind ordering; the interior FaceFlux values ARE its
off-diagonal couplings. Reflective BC adds **back-edges** `outflow(A,Ω)→inflow(B,Ω')`,
breaking strict triangularity into `L_oct + R`. The G-S vs Jacobi distinction is the
**triangular split of R**:
- **Jacobi** (current bare sweep + external `−B`): ALL reflective back-edges lagged.
  The fixed-point operator is `L⁻¹(q + Rψ_old)`. The current code is pure-Jacobi BY
  CONSTRUCTION — and does not NAME it as a split decision; it presents "reflective
  coupling is external" as an architectural fact rather than the Jacobi end of a
  spectrum.
- **Gauss-Seidel**: reflective back-edges from an already-swept face to a not-yet-swept
  face fold INTO the inverted factor. The cleanest expression: build an **(octant×face)
  reflective graph**; edges respecting the octant ordering (producer swept before
  consumer) are GS-ELIGIBLE; cycle-forming edges (true cycles) MUST be lagged. The split
  is `{e : topo-order(producer) < topo-order(consumer)}` — ONE graph computation, not a
  per-BC special case, falling directly out of the DAG machinery the sweep already uses.

Elegance payoff: structure-exposing (names interior FaceFlux as `L_oct` off-diagonal,
reflective coupling as back-edge set R; G-S/Jacobi as a graph property), algorithmic-
advantage (folding ordering-respecting reflective edges into the factor drops the
spectral radius → fewer outer iterations), structurally-simpler (one topological
comparison replaces per-geometry/per-BC branching).

First test: on a 2-D fully-reflective box, classify each reflective edge as
ordering-respecting or cycle-forming; run lag-all (Jacobi) vs fold-respecting (GS);
assert GS converges in strictly fewer outer iterations to the same fixed point.

## VERDICT on the proposed taxonomy
- "one face-flux field, boundary = edge subset, absorption = identity" IS the native
  decomposition — it is the discrete trace operator `ι*` (Frame 1).
- **Per-face standalone object: REJECTED.** Native unit = field-on-FaceLayout + zero-copy
  views, forced by THREE independent arguments: (1) vectorization (a per-face Python
  object reintroduces the 10–20× regression; the cochain frame is indifferent to storage
  granularity, so vectorization settles it for the flat buffer); (2) consistency with the
  boundary field (already flat-buffer+FaceLayout+face_view — if FaceFlux is the superset
  and the boundary is its restriction, it MUST share the pattern or ι* is not a zero-copy
  view); (3) biproduct consistency (`C¹=C¹_int⊕C¹_∂` is a biproduct; `ι*`=projection π_∂,
  `ι_*`=injection ι_∂, with `π_∂ι_∂=id` ("absorption=identity") and `π_int ι_∂=0` — holds
  cleanly only if both summands share the flat-buffer substrate so projections are slices).
- **Why interior FaceFlux is correctly ephemeral + single-role**: it is the transient
  off-diagonal of a triangular factor re-formed each sweep; the role grid
  (flux/source/residual) is a 0-cochain (CELL) concept. The boundary trace carries the
  role grid only because the boundary 1-chain is where the BC residual lives.

## CROSS-METHOD POLLINATION (SN)
- **Parallel-SN sweep scheduling (Adams/Larsen, KBA, Pautz):** the (octant×face)
  reflective graph (Frame 2) is the data structure KBA "diagonal pipelining" and Pautz
  optimal-scheduling operate on; the topological levels `i+j=k` ARE the KBA wavefront
  diagonals. The borrowing is the reflective-edge scheduling layer: schedule octants so
  ordering-respecting reflective edges are maximized (Frame 2's GS payoff). Literature
  home for the G-S/Jacobi split decision.
- **FEEC / Whitney forms:** the trace operator `ι*` is the Whitney-form trace/pullback;
  FEEC names the well-posedness of the trace (admissible BCs on a 1-cochain). Low priority
  — conceptual naming, not a code change; ORPHEUS DD is already the correct lowest-order.
- **CP / boundary integral (scattering / transfer matrix):** the reflective coupling R is
  a boundary-to-boundary scattering operator S; `T=(I−S·L⁻¹_∂)⁻¹` names the fully-reflected
  problem as a resolvent (cf. [[variant_alpha_2surface_bie_frame]]). Connects the GS split
  to the spectral radius of S.

## UNEXPLORED (confirmed low-signal for SN-sweep storage-typing — consistent set)
de Rham cohomology proper (cochain complex present but no `∂²=0` payoff — single-level DD,
B∘B≠0); MPS/tensor-networks (per-face cochain has bond dim 1 in the chain direction,
degenerate — the non-trivial MPO bond is the rank-N µ-resolved kernel,
[[trajectory_resolvent_foreign_frames]], a DIFFERENT object); symplectic (no Hamiltonian
trigger); homogenization (no scale-separation trigger); spectral-on-typing + group-theory
(re-enter at the R/scattering layer, not the storage-locus question).

Cross-refs: [[issue_168_phase_c_sweep_frame]] (faces primary, cells derived; trace law
owns the boundary DAG edge), [[issue_208_operator_algebra_frames]] (`V_bulk⊕V_boundary`
biproduct = `C¹_int⊕C¹_∂` cochain decomposition; `−B` resolvent),
[[issue_196_phase_g_chain_dag_scan_frame_attack]] (`L_oct` bidiagonal, 2×2 lower-triangular
monoid, KBA), [[variant_alpha_2surface_bie_frame]] (boundary-to-boundary scattering S).
