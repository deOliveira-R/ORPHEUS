# D5b-S2 convention crosswalk (Pattern 7) — the d≥2 UBLD wiring

The carve crosses subsystem boundaries: scheme-trait ↔ walk-storage ↔
kernel ↔ source-carrier. Written BEFORE code (coding-elegance Pattern 7).

## The `n_spatial_moments` lattice

Class-level **per-axis** 1-D basis size: `spatial_basis_per_axis: ClassVar[int]`.
- DD / Step = `1` (cell-average only).
- LD = `2` (Legendre `{1, P₁}` per axis).

Derived (d known at the call site from `len(s_axes)`):
- per-cell moments = `per_axis ** d`  (DD: 1; LD-d2: 4; LD-d3: 8).
- per-face moments = `per_axis ** (d-1)` (DD: 1; LD-d2: 2; LD-d3: 4).
- `is_multi_moment` (boolean) = `per_axis > 1` — the GATE on every wiring site.

DD/Step at `per_axis==1` ⇒ per-face=1, per-cell=1 ⇒ the rank-r face / rank-3
psi_avg stay EXACTLY (no trailing length-1 axis appended). This is the
backward-compat invariant gating every sub-step.

## Moment ordering (the load-bearing convention — from _ubld primitive)

`assemble_ubld([hx, hy], [mu_x, mu_y], ...)` uses Kronecker x-outer, y-inner:
- per-cell d=2 moment vector index = `2*o_x + o_y`, order `[bar, ŷ, x̂, x̂y]`
  (confirmed by `test_d2_exact_on_bilinear`: `psi_exact=[pbar,phat_y,phat_x,phat_xy]`).
- per-axis upstream FACE moments (transverse object), `2^{d-1}` long:
  - axis 0 (x) inflow: transverse = y → `[y-bar, y-slope]`.
  - axis 1 (y) inflow: transverse = x → `[x-bar, x-slope]`.

## Per-axis OUTGOING face (the kernel must produce, the cochain must carry)

Downstream face on axis `a` = trace of the `2^d` solution at the downstream
node on axis `a`, as a `2^{d-1}` transverse moment object.
- For the trailing axis convention used by `assemble_inflow_axis`, the
  upwind neighbour's face object on axis `a` IS the next cell's
  `upstream_face_moments` on axis `a` — so the OUT face on axis `a` must be
  in the SAME `2^{d-1}` transverse-Kronecker order the inflow on axis `a`
  consumes. This closure consistency (out-of-cell == in-of-next-cell) is
  what D5b.4 (matvec twin) + D5b.1 (round-trip) verify.

Outgoing-face reconstruction: the trace of `ψ̃(s_a=+1, transverse)` projected
onto the transverse Legendre moments. In tensor-Legendre coordinates, the
trace at node `+1` on axis `a` (P₀(+1)=1, P₁(+1)=1) sums the `o_a=0` and
`o_a=1` moment blocks: out_face_transverse_moment[t] = ψ[o_a=0, t] + ψ[o_a=1, t].

## Q_cells carrier

Today `Q_cells` reaches the kernel as `(N_oct or 1, ng, n_diag)` (scalar
average source). For LD d≥2 the kernel needs a trailing `2^d` moment axis:
`(N_oct or 1, ng, n_diag, 2^d)`. For the FLAT case (S2 — external/MMS scalar
source), only the average-moment slot (`[..., 0]`) is nonzero. The kernel
accepts BOTH:
- rank-r `Q_cells` (DD/Step/LD-scalar-source) → interpret as the average moment.
- rank-(r+1) `Q_cells` with trailing `2^d` (LD MMS moment source) → full vector.

For S2 the production sweep threads a SCALAR average source; the kernel
internally lifts it to `(…, 2^d)` with the average moment in slot 0 and the
rest zero. So NO `_SolveOperands.Q` carrier change is needed for the foundation
gates — the kernel does the lift. (A genuine MMS moment source carrier is S4.)

`_MATVEC_ZERO_SOURCE` (`np.zeros((1,1,1))`) broadcasts: the kernel reads
`Q_cells[..., None]`-style average slot, so a zero scalar lifts to zero moments.

## Domain-edge inflow (OUT-of-scope boundary — handled minimally)

`inflow[a]` per axis comes from `boundary_flux.face_view(face)[oct_idx]`,
scalar-per-face shape today. For S2 vacuum/zero domain-inflow, the cochain
seeding (`_octant_face_cochain` / `_MovingFrontier.seed`) must broadcast the
scalar inflow into a `2^{d-1}` moment-shaped ZEROS array when `is_multi_moment`.
This is the ONLY domain-edge touch (append a zeros moment axis at seeding).
If any gate forces non-zero domain inflow → STOP (S4 widens `BoundaryFlux`).

## Sites (file:line)

1. scheme.py — `DiscretizationSchemeBase`: `spatial_basis_per_axis` ClassVar + Protocol.
2. linear_discontinuous.py — `spatial_basis_per_axis=2`; `_kernel_terms` d≥2 arm
   via `_ubld`; `cell_kernel_batch`/`residual_kernel_batch` return `(psi_avg_moments,(out0,out1,...))`.
3. sweep_graph.py — `_MovingFrontier._win` alloc + `_CellSolve.cell`/`_CellResidual.cell`
   emit gated on multi-moment.
4. loss_representation.py — `_octant_face_cochain` + `_edge_outflow` moment-axis + seed zero-pad.
5. test_linear_discontinuous.py — INVERT D5b.0.
6. test_ld_ubld_primitive.py — d=2 A==A entry-wise pin (carry-forward).
