---
name: sn-sig-t-layout-drift-indexerror
description: ERR-055 — 6 SN curvilinear sweep-regression tests crashed with IndexError because they fed _sweep_1d_unified bare (nx,ng,ny) arrays after the production contract flipped to PR-INDEX-5 (ng,nx,ny); fix is test migration onto transport_sweep + PerOrdinateSource producer.
metadata:
  type: project
---

ERR-055 (2026-06-01, branch refactor/field-role-typing). Six SN
curvilinear sweep-regression tests crashed with `IndexError: index N
out of bounds for axis 1 with size 1` at `sweep_cache.py:431`
(`sig_t[:, geom.chain_idx].transpose(1,0,2)`).

**Root cause — test-side convention drift (failure mode #6), NOT a
solver bug.** The 6 tests call internal `_sweep_1d_unified(Q, sig_t,
...)` directly with bare ndarrays in the OBSOLETE `(nx, ng, ny)`
layout (e.g. `np.full((10,1,1), 0.5)` for a 10-cell ng=1 mesh). The
production contract (`sweep.py:167-168` / `:164`) is `sig_t (ng, nx,
ny)`, source `(N, ng, nx, ny)`. Introduced by commit `6cfdfd4`
(Issue #196 PR-INDEX-2) which flipped `CollisionCache` to `(ng, nx)`;
production producers (`solver.py`, `operator.py`,
`material_xs_field.py:372` `xs.sig_t.T.reshape(ng,nx,ny)`) were
migrated but these six direct-call tests were not. With ng=1, the
`(nx,1,1)` array sliced `[:,:,0]` → `(10,1)` reads as `(ng=10,
nx=1)`; chain_idx indexes 10 cells into size-1 axis → IndexError.

**Production-vs-legacy verdict:** `_sweep_1d_unified` is PRODUCTION
(sole 1-D sweep body, reached via `transport_sweep` from
`solver.py:970` + `operator.py:1946`). NOT a dying path → fix is a
test migration, not helper retirement. The green matvec twin
`test_unified_matvec_cylinder.py` already used the correct
`(ng,n_cells,1)` layout, masking the drift.

**Fix.** Migrated all 6 onto the production producer:
`PerOrdinateSource.from_isotropic(Q_iso(ng,nx,ny), sn_mesh)` swept
through `transport_sweep`. Output index swaps: scalar `phi[:,0,0] →
phi[0,:,0]`; angular `ang[n,:,0,0] → ang[n,0,:,0]`. 6/6 pass; full
test_spherical + test_cylindrical 51/51 pass; regression chunk
(unified_matvec_cylinder + native_matvec + phase_c_gates) 60 passed
6 xfailed (pre-existing).

**Methodology lessons:**
- Deterministic IndexError ⇒ skip the 8-step isolation cascade; the
  traceback names file:line directly. Token-adjacency/convention
  screen (Step 4.5) is the relevant lens.
- `coding-elegance` Pattern 7 applies to TEST fixtures: bypassing the
  production producer to hand-build bare arrays re-applies a layout
  convention at the consumer — a drift landmine. Route tests through
  the SAME producer the solver uses.
- A layout migration is incomplete until direct-`_helper`-call tests
  are migrated too (PR-INDEX migrated producers, missed 6 tests).
- ng=1 degeneracy aliases `(nx,ng,ny)` against `(ng,nx,ny)`; the swap
  only declares at the chain-index slice. ≥2 groups surfaces it
  sooner (H1 at the data-layout level).

Diagnostics deleted post-fix (covered by migrated regression gates):
`diag_01_characterize.py`, `diag_02_minimal_reproducer.py`.
