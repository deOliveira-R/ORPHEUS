---
name: dg-boundary-flux-consumer-audit
description: Pre-D-G audit of every consumer of orpheus/sn/boundary_flux.py:BoundaryFlux — production + tests, with WRITE sites (mutability) called out explicitly. Source of truth for D-G migration scope (BoundaryFlux to pure Field + SweepScratch carve).
metadata:
  type: project
---

# D-G BoundaryFlux Consumer Audit (pre-implementation)

**Date:** 2026-05-27. Plan: `.claude/plans/depth_b_field_on_function_space.md` §3.4, §6 step D-G, §11.1 invariants 5–7. Pre-D-G file: `orpheus/sn/boundary_flux.py` (301 LOC, mutable dataclass; classmethod `zeros`; views `xmin/xmax/ymin/ymax` over flat buffers; per-class dunders).

## Scope summary

- **9 production files** import or touch BoundaryFlux fields.
- **12 test files** likewise. Plus the 2-D snapshot generator.
- **Two storage shapes** are exposed today: 1-D (`xmin_face` / `xmax_face` are owning `(N, ng)` ndarrays) vs 2-D (`xmin_xmax_buf` / `ymin_ymax_buf` are owning `(N, ng, nx+1, ny)` / `(N, ng, nx, ny+1)` ndarrays carrying both BC faces AND interior wavefront cache).
- **The mutability is intentional** at every write site below (see plan §3.4 for the inversion).

## 1. Production consumers (`orpheus/...`)

### 1.1 `orpheus/sn/boundary_flux.py` (THE class — 301 LOC)
- Pre-D-G interface as described in the prompt; no behavioural change needed beyond replacing with a `Field` subclass and pushing geometry-conditional zero into `zeros_for_sn_mesh`.
- **WRITE sites inside the class**: none — all mutation comes from callers. The class itself supports mutation via field assignment (mutable `@dataclass`).
- **Migration bucket**: medium (rewrite ~ 200 LOC; Field-inherited algebra subtracts the existing dunders).

### 1.2 `orpheus/sn/sweep.py` (863 LOC, 33 refs — **HIGHEST-RISK site**)

Three sweep-body sections, each with a different write pattern.

- **1-D slab BC seed** (lines 406–426). Lazy-allocates `xmin_face` / `xmax_face` if `None` (lines 413–416) — WRITE. Then reads `bc_left_face = boundary_flux.xmin_face`, `bc_right_face = boundary_flux.xmax_face` — read.
- **1-D curvilinear BC seed** (lines 428–436). Lazy-allocates `xmax_face` if `None` (line 433) — WRITE. Then reads `bc_outer = boundary_flux.xmax_face`.
- **2-D wavefront, lines 762–858 — THE buffer split target.**
  - Lines 769–772: lazy-allocate the persistent buffers if `None` — WRITE.
  - Lines 773–774: `psi_x = boundary_flux.xmin_xmax_buf` (view, shape `(N, ng, nx+1, ny)`); `psi_y = boundary_flux.ymin_ymax_buf` (view, shape `(N, ng, nx, ny+1)`).
  - Lines 820–831: BC apply at face slot `[..., 0, :]` / `[..., nx, :]` / `[..., :, 0]` / `[..., :, ny]` — WRITE through the view into the underlying buf.
  - Lines 834–835: `psi_x_oct = psi_x[oct_idx].copy()` / `psi_y_oct = psi_y[oct_idx].copy()` — read (snapshot of incoming face cells AND interior cache for the octant).
  - Lines 840–851: `sweep_graph.apply(...)` mutates `psi_x_oct` / `psi_y_oct` in place — the per-octant arrays carry both the face flux AND the interior wavefront cache. This IS the "interior wavefront cache" that needs to migrate to SweepScratch.
  - **Lines 854–855**: `psi_x[oct_idx] = psi_x_oct`; `psi_y[oct_idx] = psi_y_oct` — **the write-back into the persistent BoundaryFlux buffer**. This is the SweepScratch migration target: only the FACE slice (`[..., 0, :]` / `[..., nx, :]` / `[..., :, 0]` / `[..., :, ny]`) is genuinely boundary state; the interior cells `[..., 1:nx, :]` and `[..., :, 1:ny]` are wavefront cache.
- **Migration bucket**: large (>100 LOC across the three sections + the buffer-split). The 2-D body is the single highest-risk section: separating face vs interior in the write-back without changing the numerics is the load-bearing carve.

### 1.3 `orpheus/sn/operator.py` (2697 LOC, 67 refs)

Multiple operator-side write sites — all of them functional-style (build a fresh `BoundaryFlux` then assign into it), well-suited to converting to immutable constructors via `replace`.

- **`transport_operator_matvec_unified` (line ~995–1090)**: reads `psi.boundary.xmax_face` / `xmin_face` (lines 1067–1068); BUILDS a fresh `m_boundary = BoundaryFlux(mesh=sn_mesh)` then writes `m_boundary.xmax_face = np.zeros((N, ng))` (line 1302), conditionally `m_boundary.xmin_face = np.zeros((N, ng))` (line 1304), then face-slice writes at lines 1308–1311 and 1316–1319 — WRITE. Functional pattern despite the mutation; D-G converts cleanly to `BoundaryFlux(values=..., layout=..., mesh=...)` constructor calls.
- **`SNStreamingOperator._apply_packed` (line 1591–1596)** — legacy 1-D path. Builds a synthetic `BoundaryFlux.zeros(sn_mesh)` then writes `boundary_legacy.xmax_face[:, :] = psi_view[:, :, -1, 0]` and conditionally `xmin_face`. **WRITE — synthetic cell-centre proxy seed**. Note: SNStreamingOperator is on retirement per plan §3.7, but the 2-D code path in this matvec body keeps the pattern alive for D-G.
- **`StreamingOperator._apply_packed_legacy` (line ~2044–2068)** — same pattern; builds `boundary_in = BoundaryFlux(mesh=sn_mesh)`; writes `xmax_face = np.zeros(...)`, conditionally `xmin_face = np.zeros(...)`, then face-slot scatter at lines 2062–2068. WRITE.
- **`StreamingOperator._apply_packed_legacy` read-back (lines 2075–2082)**: reads `result.boundary.xmax_face[eq_map.face_outer_ordinate, :]` / `xmin_face[...]` to repack — read only.
- **`StreamingOperator._apply_typed` (line ~2151)**: returns `AngularFlux(cell_values, sn_mesh, boundary=result.boundary)` — uses result.boundary as-is, no further mutation.
- **`InvertibleOperator.solve` (line ~2562–2677)**: `boundary_buf = sn_mesh.zeros_boundary_flux()` then `_copy_boundary_face_state(initial_guess.boundary, boundary_buf)` or `_copy_boundary_face_state(rhs.boundary, boundary_buf)` (lines 2660 / 2662) — WRITE via the helper below. Then passes `boundary_buf` to `transport_sweep` — sweep mutates further.
- **`_copy_boundary_face_state` helper (lines 2680–2697)**: the elementwise copy primitive — WRITES `dst.xmin_face[...] = src.xmin_face`, `dst.xmax_face[...] = src.xmax_face`, plus the two 2-D buffers. **All four write sites.** This entire helper retires in D-G (immutable BoundaryFlux → use construction-time argument passing instead of in-place buffer copies).
- **Algebra usage**: none in operator.py — no `bf + bf2` style algebra; the algebra goes through `AngularFlux + AngularFlux` which propagates via the inherited propagation chain (D-G keeps that working via Field-inherited dunders).
- **Migration bucket**: large (~300 LOC across the matvec helpers + InvertibleOperator.solve + `_copy_boundary_face_state`). Most write sites are constructor-shaped; rewire to `BoundaryFlux(...)` constructor with the values set via slice-view writes into a pre-allocated buffer OR via per-face concatenation.

### 1.4 `orpheus/sn/angular_flux.py` (513 LOC, 26 refs)

- **`__post_init__` line 134–138**: auto-allocates `BoundaryFlux.zeros(self.mesh)` when `boundary is None` — WRITE (via `object.__setattr__`). Stays in D-G but the construction call changes to `BoundaryFlux.zeros_for_sn_mesh(self.mesh)`.
- **`from_flat_with_traces` (lines 377–406)**: `boundary = BoundaryFlux.zeros(mesh)`; then conditionally writes `boundary.xmax_face[eq_map.face_outer_ordinate, :] = psi_face_outer` (line 403) and `xmin_face[...] = psi_face_inner` (line 405) — WRITE.
- **`to_flat_with_traces` (lines 408–447)**: reads `self.boundary.xmax_face[eq_map.face_outer_ordinate, :]` (line 434) and `xmin_face[...]` (line 440) — read.
- **`copy` method (lines 470–491)**: reads all four face attributes (`xmin_face`, `xmax_face`, `xmin_xmax_buf`, `ymin_ymax_buf`); constructs a fresh `BoundaryFlux(mesh=..., xmin_face=..., xmax_face=..., xmin_xmax_buf=..., ymin_ymax_buf=...)` with `.copy()` per ndarray (lines 475–485). Read-only with construction; trivial migration.
- **Algebra**: propagates to `BoundaryFlux` dunders via the `(values, boundary)` pair under the `__add__`/`__sub__`/`__mul__` chain in `angular_flux.py` (existing pattern; works unchanged with Field-inherited dunders).
- **Migration bucket**: medium (~50 LOC; `from_flat_with_traces` write block needs to become a constructor pattern with face slot scatter assembled before construction). D-H retires this file entirely, so the bucket is "as little as needed to keep tests green until D-H."

### 1.5 `orpheus/sn/solver.py` (1419 LOC, 13 refs)

- **`__init__` line 206**: `self._boundary_flux = sn_mesh.zeros_boundary_flux()` — constructor (not a write site per se; just the allocation).
- **`_solve_source_iteration` line 539**: `_copy_boundary_face_state(self._boundary_flux, q_ext_typed.boundary)` — WRITE via the helper.
- **`_solve_krylov`**: similar pattern (line 687 area; no boundary copy here — Krylov takes initial_guess via AngularFlux).
- **Migration bucket**: small (~10 LOC); `_copy_boundary_face_state` retires, replaced by passing the BoundaryFlux as an immutable construction argument to AngularFlux.

### 1.6 `orpheus/sn/solution.py` (367 LOC, 13 refs)

- All references are docstrings + the `boundary_flux` property (line 219) which delegates to `self.angular_flux.boundary` — read-only.
- **Migration bucket**: trivial (<5 LOC; docstring refresh only).

### 1.7 `orpheus/sn/geometry.py` (8 refs)

- `SNMesh.zeros_boundary_flux()` factory (lines 795–805) — delegates to `BoundaryFlux.zeros(self)`. D-G adds the new `boundary_face_layout` property used by `BoundaryFlux.zeros_for_sn_mesh(mesh)`.
- **Migration bucket**: small (~30 LOC to add `boundary_face_layout` property; replace `BoundaryFlux.zeros(self)` call).

### 1.8 `orpheus/sn/scattering.py` (1 ref)
- Single docstring mention at line 947. Trivial.

### 1.9 `orpheus/sn/axis.py` (3 refs)
- Docstrings only at lines 74, 288. Trivial.

### 1.10 `orpheus/transport/__init__.py` + `orpheus/transport/fields/__init__.py`
- Docstring placeholders for the post-D-G BoundaryFlux re-export. Trivial.

## 2. Test consumers (`tests/...`)

| File | LOC bucket | Refs | Write sites? | Notes |
|---|---|---|---|---|
| `tests/sn/test_boundary_flux_arithmetic.py` | small | 80+ | YES (10+) — `bf.xmin_face[...] = 1.0` style direct mutation throughout | The full algebra test file. Every test seeds `BoundaryFlux.zeros(sn)` then writes face buffers, then checks dunder output. Pure assignment-via-mutation pattern. D-G must rewrite to construct-then-compare or construct-fresh-with-values. ~50 LOC change. |
| `tests/sn/test_typed_fields.py` | small | 20+ | YES — `bf.xmin_xmax_buf[:, :, 0, :] = 1.0` and `bf.xmin_face[:, :] = 5.0` shape-access tests | Tests shape correctness + face slice view semantics. The 2-D `xmin_xmax_buf` writes (lines 256–271) become writes into a SweepScratch-equivalent fixture; the 1-D `xmin_face` / `xmax_face` writes (lines 284–285) become field constructor calls. ~30 LOC change. |
| `tests/sn/test_angular_flux_with_boundary.py` | medium | 30+ | YES — `psi1.boundary.xmin_face[...] = 0.5` propagation tests | Pins AngularFlux + BoundaryFlux algebra propagation. The mutations are the test fixture pattern. Re-tooling to build BoundaryFlux via constructor + replace pattern. ~40 LOC. |
| `tests/sn/test_native_matvec.py` | medium | 50+ | YES — `bf.xmax_face = np.full((N, ng), value)` etc. | Builds synthetic BoundaryFlux instances with bare construction then field-set assignment. D-G migrates to constructor with all values passed in at init. ~60 LOC. |
| `tests/sn/test_invertible_operator.py` | small | 3 | NO writes — only reads `composite_out.boundary.xmax_face` for comparison | ~5 LOC change. |
| `tests/sn/test_operators_apply_typed.py` | small | 25+ | YES — `psi.boundary.xmax_face[...] = rng.standard_normal(...)` | Random fixture pattern. ~30 LOC change. |
| `tests/sn/test_streaming_operator_decomposition.py` | small | 10+ | YES — `boundary_in.xmax_face = np.zeros(...)`, slot scatter | Same pattern as operator.py `_apply_packed_legacy`. ~15 LOC. |
| `tests/sn/test_fixed_source_g1.py` | small | 5 | NO — read-only assertions about `psi_typed.boundary.xmax_face.shape` | ~5 LOC. |
| `tests/sn/test_solution.py` | trivial | 3 | NO writes — identity assertions only | ~3 LOC. |
| `tests/sn/_test_helpers.py` | small | 10+ | YES — `bf.xmax_face = psi_view[:, :, -1, 0].copy()`; `make_boundary_flux_zero` factory | The helper `_make_boundary_with_cell_centre_proxy_seed` builds + mutates synthetic BoundaryFlux. ~20 LOC. |
| `tests/sn/test_2d_octant_sweep_equivalence.py` | medium | 15+ | YES — passes `boundary_flux` to sweep AND reads `inputs.boundary_flux.xmin_xmax_buf` post-sweep | Critical 2-D bit-identity gate. The post-sweep reads of `xmin_xmax_buf` and `ymin_ymax_buf` (lines 860, 864) move to SweepScratch reads. **~50 LOC change.** Carries `aniso_source=` argument that looks like stale legacy API — separate cleanup. |
| `tests/sn/regression/_generate_2d_octant_snapshots.py` | small | 5 | YES — reads `inputs.boundary_flux.xmin_xmax_buf` (line 85), `.ymin_ymax_buf` (line 86) after sweep | Bit-identity snapshot generator. ~10 LOC. |
| `tests/sn/test_angular_flux_with_boundary.py` | (above) | | | (counted above) |
| `tests/numerics/test_iteration.py` | trivial | 2 | NO — passes `boundary_flux = sn_mesh.zeros_boundary_flux()` to sweep | ~3 LOC. |
| `tests/numerics/test_iteration_angular_flux.py` | small | 10+ | YES — `psi.boundary.xmax_face[eq_map.face_outer_ordinate, :] = ...` | ~15 LOC. |

## 3. Sweep write-through pattern locations (the SweepScratch migration target)

The persistent 2-D buffers in `orpheus/sn/sweep.py` carry BOTH:
- **Boundary face cells**: `psi_x[:, :, 0, :]` (xmin), `psi_x[:, :, nx, :]` (xmax), `psi_y[:, :, :, 0]` (ymin), `psi_y[:, :, :, ny]` (ymax) — used by reflective-BC partners across sweep calls. Genuinely boundary trace state.
- **Interior face cells**: `psi_x[:, :, 1:nx, :]` and `psi_y[:, :, :, 1:ny]` — the wavefront cache. Conceptually scratch space the wavefront uses to thread upstream face fluxes from one cell to the next within a single sweep. Persisting these across sweep calls is incidental (the next sweep call overwrites them by upstream propagation), but the post-sweep snapshot test (`test_2d_octant_sweep_equivalence.py` lines 860, 864 and `regression/_generate_2d_octant_snapshots.py` lines 85–86) DOES read them — the bit-identity gate currently treats them as load-bearing.

**The split target is precisely the `psi_x[oct_idx] = psi_x_oct` / `psi_y[oct_idx] = psi_y_oct` write-back at sweep.py lines 854–855.** Post-D-G:
- The face-slot writes (BC apply at lines 820–831) target `BoundaryFlux.faces[face_name].values` via slice views over the flat layout — this stays in BoundaryFlux.
- The interior writes target a sweep-private `SweepScratch.psi_x[oct_idx, :, 1:nx, :]` (and similar for psi_y) — these live on SweepScratch.
- Within a single sweep call, `sweep_graph.apply` continues to consume the unified per-octant buffers (faces + interior together) — the SweepScratch is the producer; the BoundaryFlux face slots are read at apply-time and written back at end-of-sweep. Bit-identity test re-baselining: `test_2d_octant_sweep_equivalence.py` and the snapshot generator stop reading interior cells off BoundaryFlux; the snapshots either move to SweepScratch or get retired if not load-bearing for the BC contract.

Additionally, the 1-D lazy-init at lines 413–416 and 433–434 (`if xmin_face is None: ... np.zeros((N, ng))`) is a different kind of write: it's the "construct on first use" hatch the pre-D-G mutable API supports. Post-D-G immutability requires this construction logic to live in `BoundaryFlux.zeros_for_sn_mesh(mesh)` (already shipped per the plan's classmethod design) — the sweep no longer initializes; it asserts the caller passed a properly-zeroed instance.

## 4. AngularFlux.boundary access pattern (for D-H.1 sizing)

This is a D-H.1 concern (boundary moves OFF AngularFlux onto TransportState), not D-G — but the inventory feeds D-H sizing.

Files reading `.boundary.{xmax_face|xmin_face|xmin_xmax_buf|ymin_ymax_buf|xmin|xmax|ymin|ymax}` on an AngularFlux instance (9 files total):

| File | Read sites (approx) | Notes |
|---|---|---|
| `orpheus/sn/operator.py` | ~15 | matvec helpers + InvertibleOperator.solve |
| `orpheus/sn/angular_flux.py` | ~8 | self.boundary access in copy / from_flat / to_flat |
| `tests/sn/test_native_matvec.py` | ~20 | random + linearity tests |
| `tests/sn/test_operators_apply_typed.py` | ~20 | per-operator assertion sweep |
| `tests/sn/test_invertible_operator.py` | ~5 | composite operator parity |
| `tests/sn/test_streaming_operator_decomposition.py` | ~10 | decomposition equivalence |
| `tests/sn/test_fixed_source_g1.py` | ~5 | shape assertions only |
| `tests/sn/test_angular_flux_with_boundary.py` | ~15 | propagation algebra |
| `tests/numerics/test_iteration_angular_flux.py` | ~10 | iteration-protocol tests |

Estimated D-H.1 sizing: ~100 sites; rewire each from `psi.boundary.xxxx` to `state.boundary.xxxx` (where `state` is the new TransportState). Roughly 200 LOC of mechanical rename + ~50 LOC of plumbing through containers.

## 5. Estimated LOC per migration site (D-G only)

| File | Bucket | LOC est |
|---|---|---|
| `orpheus/sn/boundary_flux.py` | medium | ~200 (rewrite as Field subclass + FaceLayout) |
| `orpheus/sn/sweep.py` | large | ~120 (SweepScratch carve + 1-D lazy-init replacement + 2-D writeback split) |
| `orpheus/sn/operator.py` | large | ~150 (constructor-pattern conversion at 4 write sites + `_copy_boundary_face_state` retirement) |
| `orpheus/sn/angular_flux.py` | medium | ~50 (from_flat_with_traces / __post_init__ rewire; D-H retires this entirely) |
| `orpheus/sn/solver.py` | small | ~10 (`_copy_boundary_face_state` consumer site) |
| `orpheus/sn/solution.py` | trivial | ~5 (docstring) |
| `orpheus/sn/geometry.py` | small | ~30 (add `boundary_face_layout` property; `zeros_boundary_flux` becomes `zeros_for_sn_mesh` call) |
| `orpheus/sn/scattering.py` | trivial | ~3 (docstring) |
| `orpheus/sn/axis.py` | trivial | ~3 (docstring) |
| `orpheus/transport/__init__.py` etc. | small | ~20 (new module wiring + re-export shim) |
| `tests/sn/test_boundary_flux_arithmetic.py` | small | ~50 |
| `tests/sn/test_typed_fields.py` | small | ~30 |
| `tests/sn/test_angular_flux_with_boundary.py` | medium | ~40 |
| `tests/sn/test_native_matvec.py` | medium | ~60 |
| `tests/sn/test_operators_apply_typed.py` | small | ~30 |
| `tests/sn/test_2d_octant_sweep_equivalence.py` | medium | ~50 |
| `tests/sn/test_streaming_operator_decomposition.py` | small | ~15 |
| `tests/sn/test_fixed_source_g1.py` | trivial | ~5 |
| `tests/sn/test_invertible_operator.py` | trivial | ~5 |
| `tests/sn/test_solution.py` | trivial | ~3 |
| `tests/sn/_test_helpers.py` | small | ~20 |
| `tests/sn/regression/_generate_2d_octant_snapshots.py` | small | ~10 |
| `tests/numerics/test_iteration.py` | trivial | ~3 |
| `tests/numerics/test_iteration_angular_flux.py` | small | ~15 |
| **Total** | | **~900 LOC** |

This matches Risk-3 in plan §8 ("~15 files per audit memo" — actual count is 22 (9 prod + 13 test), but the plan's ~15 was the rough order; the LOC envelope is the load-bearing number).

## 6. Single highest-risk site

**`orpheus/sn/sweep.py:_sweep_2d_wavefront` (lines 700–858).** Specifically the persistent buffer pattern:

1. **Lines 769–774** — lazy-init of the persistent buffers + view-binding to `psi_x` / `psi_y`. Post-D-G these become reads of the BoundaryFlux face slice for the BC apply + reads of a separate `SweepScratch.psi_x` / `psi_y` for the interior.
2. **Lines 819–831** — per-octant BC apply writes face cells via `psi_x[oct_idx, :, 0, :] = full_face_x[oct_idx]` (and analogues). Post-D-G these target BoundaryFlux face values (immutable construction with the new face values).
3. **Lines 834–835** — per-octant copy of the full buffer (face + interior) into `psi_x_oct` / `psi_y_oct` for sweep-graph apply.
4. **Lines 854–855** — write-back into the persistent buffer of BOTH face AND interior. **This is the load-bearing split**: face cells are NEXT-iteration BC state (BoundaryFlux); interior cells are sweep-private scratch (SweepScratch).

The risk concentrates here because:
- The sweep is the SOLE producer of the post-sweep BoundaryFlux state today. Making BoundaryFlux immutable forces the sweep to RETURN a fresh BoundaryFlux instead of mutating the passed-in one. The signature change of `transport_sweep` ripples through `InvertibleOperator.solve`, `SNSolver._boundary_flux` plumbing, and the 2-D bit-identity tests.
- The face-vs-interior separation in the writeback isn't currently encoded — the buffer indices `[..., 0, :]` (face) vs `[..., 1:nx, :]` (interior) live implicitly in the BC apply and sweep_graph apply patterns. A `FaceLayout` descriptor needs to encode this without introducing slicing bugs.
- The 2-D snapshots in `test_2d_octant_sweep_equivalence.py` and `_generate_2d_octant_snapshots.py` read both face AND interior cells off BoundaryFlux today. The post-D-G test reads either route to SweepScratch (if the snapshot is structurally about interior wavefront) or stay on BoundaryFlux (if structurally about BC trace). Either path requires re-baselining bit-identity assertions, which is the most error-prone part of D-G.

Recommended mitigation: implement SweepScratch FIRST as a parallel scratch type (without touching BoundaryFlux), verify 2-D bit-identity on a side branch via `assert_array_equal(boundary_flux.xmin_xmax_buf, np.concatenate(face_slot_views, scratch_interior, axis=...))`, THEN flip BoundaryFlux to immutable. This is a Pattern-2 dual-view step that keeps the carve auditable.

## 7. Cross-cutting findings

- **No 1-D BoundaryFlux algebra usage in production** — `bf + bf2` style only appears in `test_boundary_flux_arithmetic.py` and propagates through AngularFlux dunders in production. The Field-inherited dunders cover both.
- **`_copy_boundary_face_state` is the canonical "in-place buffer copy" helper** (operator.py:2680). It retires in D-G — every consumer (solver.py:539, operator.py:2660/2662) switches to passing the BoundaryFlux as a construction argument.
- **`from_flat_with_traces` (angular_flux.py:377) is the most complex write site** because it scatters face values into per-ordinate slots conditional on `eq_map.face_outer_ordinate` masks. Post-D-G this becomes a builder pattern: assemble per-face ndarrays, then construct the flat buffer + layout in one call. Plan §3.4's `BoundaryFlux.zeros_for_sn_mesh(mesh)` is the canonical constructor; this matvec needs a sibling `BoundaryFlux.from_face_arrays(mesh, xmax=..., xmin=...)` factory.
- **The 2-D snapshot regression generator is a bit-identity-gate**. The plan calls out new test `tests/sn/test_sweep_scratch_split.py` (plan §6 D-G verification) — that test's purpose is precisely to pin the carve.

---

**End of audit.**
