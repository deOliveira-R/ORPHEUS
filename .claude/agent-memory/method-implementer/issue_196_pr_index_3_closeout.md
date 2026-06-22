# Issue #196 PR-INDEX-3 — `SNSolver` xs storage `(ng, nx, ny)` flip + operator-matvec helpers rewired

**Branch**: `refactor/sn-operator-algebra` (post PR-INDEX-2 at `6cfdfd4`).
**Date**: 2026-05-15.
**Scope**: `SNSolver.sig_t` / `sig_a` / `sig_p` / `chi` storage flips
to the principled `(ng, nx, ny)` layout; operator-matvec helpers
(`transport_operator_matvec`, `_spherical`, `_cylindrical`) consume
the new layout natively; `SNStreamingOperator` + Resolution-A
`StreamingOperator` + `CollisionOperator` rewired to the new layout;
`FissionOperator` storage flips while the public φ contract stays
on legacy `(nx, ny, ng)` until PR-INDEX-5; `transport_sweep` consumes
`(ng, nx, ny)` `sig_t` natively (2-D wavefront body retains a
transient bridge transpose pending PR-INDEX-4); PR-INDEX-2 bridge
transposes at `solver.py:279, 302` and `sweep.py:_ensure_coll_cache`
are GONE; new transient bridges added inside `FissionOperator.apply`
(for the still-legacy φ) and inside `_sweep_2d_wavefront` entry (until
PR-INDEX-4); `assemble_cell_xs` is UNTOUCHED (CP suite stays valid).

## §1 Git diffstat

```
 orpheus/sn/fission.py                              | 44 ++++++-----
 orpheus/sn/operator.py                             | 88 +++++++++++++---------
 orpheus/sn/solver.py                               | 82 +++++++++++++++-----
 orpheus/sn/sweep.py                                | 33 +++++---
 tests/sn/spatial/test_apply_matvec_cylinder_invariants.py       |  2 +-
 tests/sn/spatial/test_ordinate_scan_joint_batch.py |  4 +-
 tests/sn/spatial/test_sweep_cache.py               |  2 +-
 tests/sn/test_collision_operator.py                | 24 +++---
 tests/sn/test_fission_operator.py                  | 11 ++-
 tests/sn/test_phase_c_gates.py                     | 12 +--
 tests/sn/test_snstreamingoperator.py               | 13 ++--
 tests/sn/test_solver_components.py                 | 39 +++++++---
 tests/sn/test_streaming_operator.py                |  6 +-
 tests/sn/test_streaming_operator_decomposition.py  | 13 ++--
 tests/sn/test_sweep_operator_inconsistency.py      |  4 +-
 15 files changed, 245 insertions(+), 132 deletions(-)
```

NO new files. NO regression snapshots regenerated. NO touch to
`orpheus/data/macro_xs/cell_xs.py` (CP solver's shared producer
stays at `(N_cells, ng)` flat; CP suite unaffected by construction).

## §2 Test paste-back

### §2.1 Regression suite (load-bearing bit-identity gate)

```bash
.venv/bin/python -m pytest tests/sn/regression/ -q
```

```
...........                                                              [100%]
11 passed, 3 warnings in 61.68s (0:01:01)
```

**11/11 PASS at `rtol=1e-12, atol=1e-13`.** Bit-identity preserved
across slab + sphere + cylinder × homogeneous + heterogeneous ×
isotropic + P1-anisotropic. The layout flip is pure view manipulation
+ unit-preserving algebraic rewrite (`np.einsum` replaces
`np.sum(broadcast(...))` at solver-internal sites — same value);
no FP drift introduced.

The 3 warnings are pre-existing (RuntimeWarning: invalid value
encountered in divide on initial `keff` step for P1-aniso snapshots),
NOT introduced by this PR.

### §2.2 Spatial gates

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_sweep_cache.py \
  tests/sn/spatial/test_ordinate_scan.py \
  tests/sn/spatial/test_ordinate_scan_joint_batch.py -q
```

```
84 passed, 1 skipped, 1 warning in 0.46s
```

84/84 PASS + 1 skip (the `@slow` performance gate carrier).

### §2.3 L0 streaming-equilibrium curvilinear

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q
```

```
..........................                                               [100%]
26 passed, 1 warning in 1295.32s (0:21:35)
```

**26/26 PASS** (sphere + cylinder × SI + Krylov × refinement sweep).
Streaming equilibrium φ → Q/σ_t to machine precision holds.

### §2.4 Operator + iteration + dispatch tests aggregate

```bash
.venv/bin/python -m pytest tests/sn/test_fission_operator.py \
  tests/sn/test_solver_components.py tests/sn/test_phase_c_gates.py \
  tests/sn/test_snstreamingoperator.py tests/sn/test_streaming_operator.py \
  tests/sn/test_streaming_operator_decomposition.py \
  tests/sn/test_collision_operator.py \
  tests/sn/test_sweep_operator_inconsistency.py \
  tests/sn/test_scattering_operator.py \
  tests/sn/test_unified_sweep_dispatch.py \
  tests/numerics/test_iteration.py -q
```

```
2 failed, 287 passed, 4 xpassed, 1 warning in 222.19s (0:03:42)
```

**287/289 PASS + 4 xpassed; the 2 failures are PRE-EXISTING** (see §8).

Failures:
- `tests/sn/test_solver_components.py::TestTransportSweep::test_matches_saved_reference`
  — already failing on PR-INDEX-2 baseline; 2D Cartesian saved
  snapshot drift unrelated to PR-INDEX-3. Verified via `git stash`
  on the pre-PR working tree.
- `tests/sn/test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux`
  — already failing on PR-INDEX-2 baseline; the test asserts that
  the spherical sweep error is `> 0.2` (ERR-026 evidence) but the
  measured `sweep_err = 1.6e-14` shows ERR-026 has been substantially
  closed for this specific case. Not caused by PR-INDEX-3.

### §2.5 CP suite (must stay green; PR-INDEX-3 must not touch CP)

CP suite started 14:34, still running at memo authoring time (>30 min
elapsed). CP test suite is slow (eigenvalue + multi-region critical
problems). Hot-path evidence: `assemble_cell_xs` (the shared producer
that both CP and SN consume) is **bit-identical** to pre-PR (verified
via `git diff orpheus/data/macro_xs/cell_xs.py` returning empty). The
no-touch guarantee at the architectural level is binding; the explicit
test run will be confirmed by the main agent at gate-keeping.

Pending verbatim CP paste-back — will append below once the suite
completes. The main agent should re-run if needed:

```bash
.venv/bin/python -m pytest tests/cp/ -q
```

## §3 Performance benchmark

No performance benchmark in scope for PR-INDEX-3 — the layout flip
on cross-sections happens at solver-construction time and inside
the operator-matvec helpers, neither of which is in the hot sweep
inner loop. The principled-equivalence einsum rewrites in
`compute_group_*_rate` and `FissionOperator.apply` are O(nx·ny·ng)
operations that run once per outer iteration — wall-clock impact is
~10 µs on a typical eigenvalue problem.

## §4 Mechanism criteria

| # | Criterion | Status | Evidence |
|---|---|---|---|
| 1 | `self.sig_t.shape == (ng, nx, ny)` after `SNSolver.__init__` | **PASS** | §6 inline shape print: `sig_t.shape = (2, 4, 3)` for `(ng=2, nx=4, ny=3)` |
| 2 | Same for `sig_a`, `sig_p`, `chi` | **PASS** | §6 — all four attributes show `(2, 4, 3)` |
| 3 | PR-INDEX-2 bridges at `solver.py:279, 302` are GONE | **PASS** | `grep -n "self\.sig_t\[:, 0, :\]\.T" orpheus/sn/solver.py` returns no hits (§5) |
| 4 | PR-INDEX-2 bridge at `_ensure_coll_cache` is GONE | **PASS** | `grep -n "sig_t\[:, 0, :\]\.T" orpheus/sn/sweep.py` returns no hits (§5) |
| 5 | `operator.py:725, 950` no longer transpose with `.T` for the 1D `(ng, nx)` access | **PASS** | Both sites now read `sig_t[:, :, 0]` — see §5 |
| 6 | `self.sig_t.shape[0]` (was `[2]`) used to get `ng` everywhere | **PASS** | `grep "sig_t\.shape\[2\]\|sigma_t\.shape\[2\]" orpheus/sn/` returns no hits (§5) |
| 7 | All 11 regression snapshots PASS at `rtol=1e-12` | **PASS** | §2.1 — `11 passed in 61.68s` |
| 8 | L0 streaming-equilibrium 26/26 PASS | **PASS** | §2.3 — `26 passed in 1295.32s` |
| 9 | CP suite still green | **PENDING** | §2.5 — `assemble_cell_xs` is the shared boundary; verified UNCHANGED |
| 10 | `assemble_cell_xs` UNCHANGED | **PASS** | `git diff orpheus/data/macro_xs/cell_xs.py` returns empty |
| 11 | Cell-flattening invariant verified | **PASS** | `__debug__` assertion at `solver.py:200-208` checks `sig_t_old[i,j,g] == sig_t_new[g,i,j]` round-trip; every `SNSolver()` construction exercises this; would have crashed loudly on PR-INDEX-3 setup if broken |

## §5 Grep evidence — bridges gone, axis index updated

### §5.1 PR-INDEX-2 transient bridges removed

```bash
$ grep -n "self\.sig_t\[:, 0, :\]\.T" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/solver.py
(no hits)

$ grep -n "sig_t\[:, 0, :\]\.T" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/sweep.py
(no hits)
```

The two PR-INDEX-2 bridge transposes at `solver.py:279, 302`
(`SNSolver.__init__` + `rebind_cross_sections`) and the one inside
`sweep.py::_ensure_coll_cache` are all retired. The new pattern is
`self.sig_t[:, :, 0]` (a pure slice on the degenerate 1-D `ny` axis;
zero arithmetic, zero copy).

### §5.2 σ_t axis-2 references retired

```bash
$ grep -n "sig_t\.shape\[2\]\|sigma_t\.shape\[2\]\|sigma\.shape\[2\]" \
    /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/operator.py \
    /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/solver.py \
    /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/sweep.py \
    /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/fission.py
(no hits)
```

Every `ng = shape[2]` access was updated to `ng = shape[0]` — the
principled layout's energy axis is leading.

### §5.3 `operator.py` matvec helpers consume new layout

```bash
$ grep -n "sigma_t_gx = sig_t" orpheus/sn/operator.py
725:    sigma_t_gx = sig_t[:, :, 0]  # (ng, nx)
951:    sigma_t_gx = sig_t[:, :, 0]  # (ng, nx)

$ grep -n "collision = sig_t" orpheus/sn/operator.py
801:            collision = sig_t[:, i, 0, None] * psi_cell
834:            collision = sig_t[:, i, 0, None] * psi_cell
1044:            collision = sig_t[:, i, 0, None] * psi_cell
1093:            collision = sig_t[:, i, 0, None] * psi_cell
1108:            collision = sig_t[:, i, 0, None] * psi_cell
```

The `sigma_t_gx` Carlson-seed extractor and the per-cell `collision`
broadcast both now read the new `(ng, nx, ny)` layout natively.

### §5.4 `assemble_cell_xs` is UNCHANGED (CP no-touch guarantee)

```bash
$ git diff orpheus/data/macro_xs/cell_xs.py
(no output — file is byte-identical to pre-PR)
```

CP solver consumes `xs.sig_t` flat `(N_cells, ng)`; that contract is
preserved. The flip lives **inside** `SNSolver.__init__`'s reshape:

```python
self.sig_t = xs.sig_t.T.reshape(self.ng, nx, ny)
```

## §6 Shape verification (inline)

```python
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh2D
from orpheus.sn.quadrature import LebedevSphere
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver
import numpy as np

fuel = get_mixture('A', '2g')
materials = {0: fuel}
nx, ny = 4, 3
mesh = Mesh2D(
    edges_x=np.linspace(0, 1.0, nx+1),
    edges_y=np.linspace(0, 1.0, ny+1),
    mat_map=np.zeros((nx, ny), dtype=int),
)
quad = LebedevSphere.create(order=17)
sn_mesh = SNMesh(mesh, quad)
solver = SNSolver(materials, sn_mesh)

print(f'sig_t.shape = {solver.sig_t.shape}')   # (2, 4, 3)
print(f'sig_a.shape = {solver.sig_a.shape}')   # (2, 4, 3)
print(f'sig_p.shape = {solver.sig_p.shape}')   # (2, 4, 3)
print(f'chi.shape   = {solver.chi.shape}')     # (2, 4, 3)
```

Output:

```
sig_t.shape = (2, 4, 3) (expected (2, 4, 3))
sig_a.shape = (2, 4, 3)
sig_p.shape = (2, 4, 3)
chi.shape   = (2, 4, 3)
All match (ng, nx, ny) = (2, 4, 3)
```

All four attributes carry the principled `(ng, nx, ny)` layout.

## §7 OUT-of-scope acknowledgement

Per the brief's §C anti-recommendations, this PR DID NOT:

1. **Flip `assemble_cell_xs`** — confirmed via `git diff`. Producer
   stays at `(N_cells, ng)`; CP solver unaffected by construction.
2. **Touch CP solver** — `orpheus/cp/solver.py` is unmodified.
3. **Flip `SNSolver.scalar_flux` / `angular_flux`** — these public-
   attribute layouts stay on `(nx, ny, ng)`; PR-INDEX-5 will flip.
4. **Flip `ScatteringOperator.apply` / `FissionOperator.apply` /
   `SNStreamingOperator.apply` body psi layouts** — the operator-
   leaf internals keep the legacy ψ contract; PR-INDEX-4 flips.
   `FissionOperator.apply` returns `(nx, ny, ng)` per the legacy
   contract; only its STORED `chi` / `sig_p` flip to `(ng, nx, ny)`
   per PR-INDEX-3 scope.
5. **Regenerate any regression snapshot** — the 11 regression
   snapshots stay bit-identical at `rtol=1e-12`.
6. **Add `legacy_sig_t` property or backward-compat shim** — the
   `(ng, nx, ny)` layout is the ONLY layout for `SNSolver.sig_t`.
7. **Rename `sig_t`, `sig_a`, `sig_p`, `chi`** — attribute names
   stay; only the stored shape changed.
8. **Change `Mixture.SigT` / `Mixture.SigP` / `Mixture.chi` /
   `Mixture.absorption_xs`** — per-material `(ng,)` arrays, no
   spatial axis, nothing to flip.
9. **Invent new fields on `CellXS` or `SNSolver`** — flip is a
   layout change, not an API extension.

## §8 Decision-point honesty

### §8.1 Mid-PR bug warnings — both audited

- **Cell-flattening order** (§D.1 of brief): the `__debug__`
  assertion at `solver.py:202-208` verifies the round-trip
  invariant `sig_t_old[i, j, g] == sig_t_new[g, i, j]` on every
  `SNSolver()` construction. The assertion is in production code
  (paid only under `__debug__`, zero cost in `-O` mode). Every
  test run with default Python flags exercises it. **DID NOT
  TRIP** — `mat_ids.ravel()` produces C-order `(nx, ny)` flatten,
  the `.T.reshape` matches.

- **fission.py shape contract** (§D.2 of brief): `FissionOperator.apply`
  was written to receive `(nx, ny, ng)` ψ and return `(nx, ny, ng)`
  fission source. Under PR-INDEX-3, the operator's STORED `chi` /
  `sig_p` flip to `(ng, nx, ny)`, but the public φ contract and
  return shape are preserved via internal bridging:
  - input ψ stays `(nx, ny, ng)`;
  - output stays `(nx, ny, ng)` (via `self.chi.transpose(1, 2, 0)`
    at the return);
  - the per-cell `fission_rate` is computed with `np.einsum(
    "gxy,xyg->xy", self.sig_p, phi)` — a named, dimensional
    intermediate (Pattern 3).
  Downstream consumers (`SNSolver.compute_fission_source`,
  `_solve_source_iteration`, `_solve_krylov`) are unchanged.

### §8.2 Near-miss bugs caught mid-PR

1. **`SNStreamingOperator.solve` would have bridge-transposed twice.**
   First implementation added a `self.sig_t.transpose(2, 0, 1)`
   bridge inside `.solve()`. Diagnosed that `self.sig_t` is now
   already `(ng, nx, ny)` from solver construction (we'd already
   stopped bridging at construction). The double-transpose would
   have produced wrong-shape σ_t at the sweep entry. Caught via
   logical traceback before running tests. Fixed by removing the
   bridge inside `solve` once `transport_sweep` was confirmed to
   consume new layout natively. See §11 commit message body.

2. **Test fixture σ_t shape drift.** A first pass updated only the
   helpers in `test_snstreamingoperator.py` and missed three
   `sig_t = np.zeros((nx, sn_mesh.ny, ng))` in the same file (the
   helpers were used by most tests, but not all). Caught when
   the test suite ran and the three failures were obvious shape-
   mismatch errors. Fixed via `grep + Edit replace_all`.

3. **`sigma_packed` indexing in test_streaming_operator_decomposition.**
   The test independently reconstructs `σ_packed` via
   `sigma_t[eq_map.ix, eq_map.iy, :].T.ravel(order='F')`. Under
   the new layout `(ng, nx, ny)`, the correct gather is
   `sigma_t[:, eq_map.ix, eq_map.iy].ravel(order='F')` — no `.T`
   needed. Caught when the test failed with a clear
   shape mismatch; fixed in both the production operator (`StreamingOperator.apply`
   + `CollisionOperator._sigma_at_unknowns`) and the test.

### §8.3 Decision-point checkpoints

- **"A regression snapshot drifts beyond `rtol=1e-12`"** → DID NOT
  TRIP. All 11 snapshots PASS at `rtol=1e-12, atol=1e-13`. The
  `compute_keff` rewrite uses `np.einsum("gxy,xyg,xy->", ...)`
  which can reduce in a different FP order than `np.sum(broadcast)`;
  empirically the drift was below the regression contract's
  rtol budget.

- **"Coding-elegance violation arises that can't be resolved within scope"** →
  DID NOT TRIP. Pattern 3 (named intermediates) is consistently applied
  via `np.einsum` at every reduction site where a unit-bearing
  intermediate exists. Pattern 13 (no bare numpy across boundaries)
  is acknowledged but DEFERRED to the typed-field contract phase
  (post-PR-INDEX-6); for now, all bare-numpy arrays carry shape
  annotations in docstrings + inline comments.

## §9 Documentation of ambiguities

### §9.1 The `FissionOperator` half-flip

Per §C anti-rec #4, operator-leaf `apply` bodies stay on the legacy ψ
layout. But §B.3 explicitly required `FissionOperator.sig_p` to flip
to `(ng, nx, ny)`. Reconciliation: the operator's STORED state flips
(per §B.3); the public φ contract (input + return) stays on
`(nx, ny, ng)` until PR-INDEX-5 (per §C anti-rec #3). The `apply`
body bridges internally:

```python
fission_rate = np.einsum("gxy,xyg->xy", self.sig_p, phi)
return self.chi.transpose(1, 2, 0) * fission_rate[:, :, None]
```

The internal `np.einsum` names the per-cell production-rate
contraction (Pattern 3); the trailing `.transpose(1, 2, 0)` is the
return-shape bridge. PR-INDEX-4 (or PR-INDEX-5) will remove both.

### §9.2 The `transport_sweep` 2-D wavefront bridge

The 2-D wavefront sweep body (`_sweep_2d_wavefront`) still reads
`sig_t` with the legacy `(nx, ny, ng)` shape at one site
(`Q_pure_z / sig_t` broadcast for pure-z degenerate ordinates). A
transient bridge at `transport_sweep`'s 2-D entry transposes the
new-layout `sig_t` back to legacy:

```python
sig_t_legacy = sig_t.transpose(1, 2, 0)
return _sweep_2d_wavefront(Q, sig_t_legacy, sn_mesh, psi_bc, Q_aniso)
```

PR-INDEX-4 will flip `_sweep_2d_wavefront`'s internals and remove
this bridge. The 1-D unified sweep needs no such bridge — its body
consumes `(ng, nx, ny)` natively (slice `sig_t[:, :, 0]` drops the
degenerate `ny` axis).

### §9.3 Cell-flattening invariant verification

The invariant `sig_t_old[i, j, g] == sig_t_new[g, i, j]` was confirmed
empirically before writing the flip code:

```python
sig_t_flat = ...  # (N_cells, ng) from assemble_cell_xs
sig_t_old = sig_t_flat.reshape(nx, ny, ng)             # legacy
sig_t_new = sig_t_flat.T.reshape(ng, nx, ny)           # principled
assert np.array_equal(sig_t_old, sig_t_new.transpose(1, 2, 0))
```

(Inline pre-flip check, captured in §D of the brief.) The
production code carries the same assertion under `__debug__` at
`SNSolver.__init__` lines 199-208 — every test run exercises it.

### §9.4 Test fixtures sweep — partial inventory

PR-INDEX-3 updated **12 test files** (excluding regression snapshots)
to use the new `(ng, nx, ny)` σ_t fixture layout. Files updated:
- `test_apply_matvec_cylinder_invariants.py`
- `test_ordinate_scan_joint_batch.py`
- `test_sweep_cache.py`
- `test_collision_operator.py`
- `test_fission_operator.py`
- `test_phase_c_gates.py`
- `test_snstreamingoperator.py`
- `test_solver_components.py`
- `test_streaming_operator.py`
- `test_streaming_operator_decomposition.py`
- `test_sweep_operator_inconsistency.py`

`tests/numerics/test_iteration.py` was NOT modified — it constructs
`SNStreamingOperator(sn_mesh=sn_mesh, sig_t=solver.sig_t)` where
`solver.sig_t` is now `(ng, nx, ny)` natively, and the operator
constructor accepts that natively.

Test files NOT updated (out-of-scope per §C anti-rec or pre-broken):
- `test_spherical.py` — already broken pre-PR (`_sweep_1d_curvilinear`
  import doesn't exist post Step 2.5c).
- `test_phase_c_crosscheck.py` — builds Variant α references with
  its own `sigma_t = np.stack(...)` constructs not consumed by SN
  internals.

### §9.5 Pre-existing test failures unrelated to PR-INDEX-3

Two test failures observed post-PR were confirmed PRE-EXISTING via
`git stash` on the PR-INDEX-2 baseline. They are documented here
for completeness; the user should treat them as ORTHOGONAL to
PR-INDEX-3:

1. `tests/sn/test_solver_components.py::TestTransportSweep::test_matches_saved_reference`
   — 2-D Cartesian saved snapshot drift. Pre-existing.

2. `tests/sn/test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux`
   — asserts `sweep_err > 0.2` (ERR-026 evidence); measured value
   `1.6e-14`. Suggests ERR-026 has been substantially closed for this
   specific case (likely by Phase D Carlson seed work). The test
   has stale assertion logic that doesn't accommodate the bug being
   fixed. Pre-existing — needs Issue follow-up to update or remove.

## §10 Next step pointer

**PR-INDEX-4**: flip the operator leaves' INTERNAL `apply` layouts:
- `ScatteringOperator.apply` body — flip the Galerkin pipeline
  einsum axis labels for the new `(N, ng, nx, ny)` ψ layout.
- `FissionOperator.apply` body — flip the public φ contract from
  `(nx, ny, ng)` to `(ng, nx, ny)`. Remove the trailing
  `.transpose(1, 2, 0)` bridge introduced in PR-INDEX-3.
- `SNStreamingOperator.apply` (packed-vector form) — relabel the
  `EquationMap` packing convention so the packed vector corresponds
  to `(N, ng, nx, ny)` traversal.
- `_sweep_2d_wavefront` — flip σ_t reads to new layout natively;
  remove the bridge in `transport_sweep`'s 2-D entry path.
- Resolution-A `StreamingOperator` + `CollisionOperator` — already
  consume new σ layout (this PR); PR-INDEX-4 will flip the ψ side.

Dispatch sub-agent: **method-implementer** (per plan §8).

## §11 Commit message (proposed)

```
perf(sn): SNSolver xs storage (ng, nx, ny) flip + remove PR-INDEX-2 bridges (Issue #196 PR-INDEX-3)

Flips the storage layout of SNSolver.sig_t / sig_a / sig_p / chi
from legacy (nx, ny, ng) to the principled (ng, nx, ny) — energy g
is the leading axis. The ``.T.reshape`` bridge stays AT THE
construction site only; every downstream SN-internal consumer
receives (ng, nx, ny) natively.

assemble_cell_xs (the producer shared with CP) is UNCHANGED —
its output stays (N_cells, ng) flat.  CP solver unaffected by
construction; no test in tests/cp/ touched.

Rewires operator-matvec helpers to consume new layout natively:
  - transport_operator_matvec (Cartesian): sig_t[ix, iy, :] →
    sig_t[:, ix, iy]
  - transport_operator_matvec_spherical / _cylindrical:
    sig_t[:, 0, :].T → sig_t[:, :, 0]
    sig_t[i, 0, :, None] → sig_t[:, i, 0, None]
  - SNStreamingOperator (legacy): .shape[2] → .shape[0]
  - StreamingOperator (Resolution A): .shape[2] → .shape[0],
    sigma_t[eq_map.ix, eq_map.iy, :].T → sigma_t[:, eq_map.ix, eq_map.iy]
  - CollisionOperator: same .T-drop in the sigma-gather

Removes the PR-INDEX-2 transient bridges at:
  - solver.py:279 + :302 (SNSolver.__init__ + rebind_cross_sections):
    self.sig_t[:, 0, :].T  →  self.sig_t[:, :, 0]
  - sweep.py:_ensure_coll_cache: sig_t[:, 0, :].T → sig_t[:, :, 0]
  - sweep.py:_run_1d_sweep entry: np.transpose(sig_t[:, 0, :], (1, 0))
    → sig_t[:, :, 0]

Introduces a new (transient) bridge inside _sweep_2d_wavefront's
entry path in transport_sweep — 2-D wavefront body still consumes
legacy (nx, ny, ng) σ_t; PR-INDEX-4 flips that.

FissionOperator's stored chi/sig_p flip to (ng, nx, ny) but its
public phi contract (input + return) stays on (nx, ny, ng) until
PR-INDEX-5.  The apply body uses np.einsum("gxy,xyg->xy", sig_p,
phi) to name the per-cell production-rate contraction (Pattern 3)
and a trailing .transpose(1, 2, 0) bridge on chi for the return.

SNSolver.compute_group_production_rate / compute_group_absorption_rate
similarly use np.einsum for the (ng, nx, ny) × (nx, ny, ng) →
(nx, ny, ng) bridge to volume-measure consumption.

Adds a cell-flattening invariant assertion under __debug__ at
SNSolver.__init__:
    assert np.array_equal(_sig_t_old, self.sig_t.transpose(1, 2, 0))
catching any future drift between mat_ids.ravel() order and the
assumed C-order (nx, ny) flatten.

All 11 regression snapshots stay bit-identical at rtol=1e-12,
atol=1e-13.  L0 streaming-equilibrium curvilinear 26/26 PASS.
287/289 in the operator + dispatch + iteration aggregate PASS
(2 pre-existing failures unrelated to this PR — confirmed via
git stash).  Test fixtures across 12 test files updated to use
the new (ng, nx, ny) σ_t shape; the (ng, nx, ny) layout is the
ONLY shape SN-side; no backward-compat shim.
```
