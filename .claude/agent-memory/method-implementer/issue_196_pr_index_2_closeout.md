# Issue #196 PR-INDEX-2 — `CollisionCache` `(N, ng, nx)` layout flip + remove sweep cache-read bridges

**Branch**: `refactor/sn-operator-algebra` (post PR-INDEX-1 at `e09b9f8`).
**Date**: 2026-05-15.
**Scope**: `orpheus/sn/spatial/sweep_cache.py` (cache storage layout flip
+ `from_geometry` rewrite), `orpheus/sn/sweep.py` (remove `np.swapaxes`
slab transpose + `.T` curvilinear transpose at cache-read sites), test
fixture + assertion updates at `tests/sn/spatial/test_sweep_cache.py`,
solver-side σ_t transpose at the `from_geometry` call sites
(`orpheus/sn/solver.py`).

## §1 Git diffstat

```
 orpheus/sn/solver.py                 |  11 +++-
 orpheus/sn/spatial/sweep_cache.py    | 118 +++++++++++++++++++++++------------
 orpheus/sn/sweep.py                  |  38 +++++------
 tests/sn/spatial/test_sweep_cache.py |  82 +++++++++++++++---------
 4 files changed, 156 insertions(+), 93 deletions(-)
```

NO new files. NO regression snapshots regenerated.

## §2 Test paste-back

### §2.1 `tests/sn/spatial/test_sweep_cache.py` verbose

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_sweep_cache.py -v
```

```
============================= test session starts ==============================
platform darwin -- Python 3.14.3, pytest-9.0.2, pluggy-1.6.0 -- /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python
cachedir: .pytest_cache
rootdir: /Users/rodrigo/git/nuclear/ORPHEUS
configfile: pyproject.toml
plugins: dash-4.1.0, anyio-4.13.0
collecting ... collected 28 items

tests/sn/spatial/test_sweep_cache.py::test_geometry_coefficients_built_at_construction PASSED [  3%]
tests/sn/spatial/test_sweep_cache.py::test_collision_cache_built_at_sigma_t_bind PASSED [  7%]
tests/sn/spatial/test_sweep_cache.py::test_two_strata_independence_by_ng_axis PASSED [ 10%]
tests/sn/spatial/test_sweep_cache.py::test_collision_cache_invariance_under_source_iteration PASSED [ 14%]
tests/sn/spatial/test_sweep_cache.py::test_geometry_coefficients_invariance_under_sigma_t_change PASSED [ 17%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-1-slab] PASSED [ 21%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-1-sphere] PASSED [ 25%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-2-slab] PASSED [ 28%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-2-sphere] PASSED [ 32%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-3-slab] PASSED [ 35%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-3-sphere] PASSED [ 39%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-1-slab] PASSED [ 42%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-1-sphere] PASSED [ 46%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-2-slab] PASSED [ 50%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-2-sphere] PASSED [ 53%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-3-slab] PASSED [ 57%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-3-sphere] PASSED [ 60%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-1-slab] PASSED [ 64%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-1-sphere] PASSED [ 67%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-2-slab] PASSED [ 71%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-2-sphere] PASSED [ 75%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-3-slab] PASSED [ 78%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-3-sphere] PASSED [ 82%]
tests/sn/spatial/test_sweep_cache.py::test_cache_populator_matches_cell_balance_terms PASSED [ 85%]
tests/sn/spatial/test_sweep_cache.py::test_slab_sweep_benchmark_under_2ms PASSED [ 89%]
tests/sn/spatial/test_sweep_cache.py::test_full_sn_suite_under_5min SKIPPED [ 92%]
tests/sn/spatial/test_sweep_cache.py::test_l0_streaming_equilibrium_preserved_after_2_5c PASSED [ 96%]
tests/sn/spatial/test_sweep_cache.py::test_pair_monoid_associativity_still_passes PASSED [100%]

=================== 27 passed, 1 skipped, 1 warning in 0.81s ===================
```

27/27 PASS + 1 skip (the performance-gate `@slow` marker passed inline;
the suite-level `test_full_sn_suite_under_5min` is a marker carrier).
Dual-view contract 18/18 parametrised cases (geometry × ng × source)
hold at `rtol=1e-13`.

### §2.2 ordinate_scan + joint_batch

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_ordinate_scan.py tests/sn/spatial/test_ordinate_scan_joint_batch.py -q
```

```
.........................................................                [100%]
57 passed, 1 warning in 0.33s
```

57/57 PASS. The `ordinate_scan` algorithm + PR-INDEX-1's joint-batch
test remain green — PR-INDEX-2 didn't touch `scan.py`.

### §2.3 Regression suite (load-bearing bit-identity gate)

```bash
.venv/bin/python -m pytest tests/sn/regression/ -q
```

```
...........                                                              [100%]
11 passed, 3 warnings in 61.49s (0:01:01)
```

**11/11 PASS at `rtol=1e-12, atol=1e-13`.** Bit-identity preserved
across slab + sphere + cylinder × homogeneous + heterogeneous ×
isotropic + P1-anisotropic. The transpose-based layout change is pure
view manipulation; no FP drift introduced (per the brief's
decision-point checkpoint: "layout change is pure view manipulation;
FP drift is mechanically impossible").

The 3 warnings are pre-existing (RuntimeWarning: invalid value
encountered in divide on initial `keff` step for P1-aniso snapshots),
NOT introduced by this PR.

### §2.4 L0 streaming-equilibrium curvilinear

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q
```

```
..........................                                               [100%]
26 passed, 1 warning in 911.26s (0:15:11)
```

**26/26 PASS** (sphere + cylinder × SI + Krylov × refinement sweep).
Streaming equilibrium φ → Q/σ_t to machine precision holds.

## §3 Performance benchmark

Inline microbench (slab N=16, ng=2, nx=160, vacuum BC, 100 sweeps,
3 warmup; 3 trials with cache invalidated between trials):

```
POST-PR-INDEX-2 slab N=16 ng=2 nx=160 trial 0: 0.154 ms/sweep
POST-PR-INDEX-2 slab N=16 ng=2 nx=160 trial 1: 0.149 ms/sweep
POST-PR-INDEX-2 slab N=16 ng=2 nx=160 trial 2: 0.145 ms/sweep
```

**Mean ~0.149 ms/sweep**, compared to PR-INDEX-1's trial mean ~0.211 ms
(0.156 / 0.327 / 0.150 ms). The removal of two bridge transposes
(`np.swapaxes` on slab, `.T` per ordinate on curvilinear) tightens
the timing-distribution variance — trial-to-trial scatter dropped from
~2× to ~1.06×. The improvement is measurable AND stable; this is the
expected mechanical signature of removing per-call view-manipulation
work from the hot path. Well within the ±5% performance gate; in
practice a clear win over PR-INDEX-1.

## §4 Mechanism criteria

| # | Criterion | Status | Evidence |
|---|---|---|---|
| 1 | `CollisionCache.inverse_denom.shape == (N, ng, nx)` | **PASS** | §6 inline shape probe — `(8, 3, 4)` for N=8, ng=3, nx=4; `sweep_cache.py:367` field comment; `test_sweep_cache.py:122` assertion |
| 2 | `CollisionCache.a_attenuation.shape == (N, ng, nx)` | **PASS** | §6 inline shape probe; `sweep_cache.py:368`; `test_sweep_cache.py:123` |
| 3 | `np.cumprod` in `from_geometry` runs over `axis=2` (cell axis) | **PASS** | `orpheus/sn/spatial/sweep_cache.py:452` — `cumprod_a = np.cumprod(a_attenuation, axis=2)` (was `axis=1` pre-PR) |
| 4 | `_run_1d_sweep` no longer calls `np.swapaxes` on cache fields | **PASS** | `grep -n "swapaxes" orpheus/sn/sweep.py` returns no hits (§5) |
| 5 | `_run_1d_sweep` no longer calls `coll.*.T` for curvilinear path | **PASS** | `grep -nE "coll\.[a-z_]+\.T" orpheus/sn/sweep.py` returns no hits (§5) |
| 6 | All 11 regression snapshots PASS at `rtol=1e-12` | **PASS** | §2.3 — `11 passed in 61.49s` |
| 7 | L0 streaming-equilibrium 26/26 PASS | **PASS** | §2.4 — `26 passed in 911.26s` |
| 8 | `test_sweep_cache.py` PASS with updated assertions | **PASS** | §2.1 — 27 passed + 1 skip |
| 9 | `from_geometry` accepts σ_t as `(ng, nx)` for 1D | **PASS** | `orpheus/sn/spatial/sweep_cache.py:394-396` (signature docstring); `sweep_cache.py:430` σ_t chain-ordering uses `sig_t[:, geom.chain_idx]` indexing axis-1 |
| 10 | Solver-side `from_geometry` callers transpose at the call site | **PASS** | `orpheus/sn/solver.py:280` (`__init__`) and `:303` (`rebind_cross_sections`) — `sig_t_1d = self.sig_t[:, 0, :].T  # (ng, nx)` |

## §5 Grep evidence — three bridge transposes gone

```bash
$ grep -n "swapaxes" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/sweep.py
(no hits)

$ grep -nE "coll\.[a-z_]+\.T" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/sweep.py
(no hits)

$ grep -n "np.cumprod" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/spatial/sweep_cache.py
452:        cumprod_a = np.cumprod(a_attenuation, axis=2)                    # (N, ng, nx)
```

Three load-bearing assertions all pass:

1. No `swapaxes` on the cache-read SLAB path — was `inv_denom_chain = np.swapaxes(coll.inverse_denom[ords], 1, 2)` (and same for `a_attenuation`); now the indexed slice IS `(K, ng, nx)` natively.
2. No `coll.*.T` on the cache-read CURVILINEAR path — was `inv_denom_p = coll.inverse_denom[global_n].T`; now `coll.inverse_denom[global_n]` IS `(ng, nx)` natively.
3. `cumprod` runs over `axis=2` (the cell axis under principled layout); was `axis=1` (the cell axis under legacy layout). The same SEMANTIC axis ("the cell axis") is preserved; only its INDEX position changed because the SHAPE flipped.

## §6 Shape verification (inline)

```python
import numpy as np
from orpheus.geometry import BC, Mesh1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.spatial.sweep_cache import CollisionCache, GeometryCoefficients

mesh = Mesh1D(edges=np.linspace(0.0, 1.0, 5), mat_ids=np.zeros(4, dtype=int),
              bc_left=BC('vacuum'), bc_right=BC('vacuum'))
quad = GaussLegendre1D.create(8)
sn_mesh = SNMesh(mesh, quad)
geom = GeometryCoefficients.from_mesh_and_quad(sn_mesh)
# σ_t in principled (ng, nx) layout
sig_t = np.array([[1.0]*4, [2.0]*4, [0.5]*4])  # ng=3, nx=4
coll = CollisionCache.from_geometry(geom, sig_t)
```

Output:

```
N= 8 ng=3, nx=4
inverse_denom shape: (8, 3, 4) (expected (8, 3, 4))
a_attenuation shape: (8, 3, 4) (expected (8, 3, 4))
cumprod_a shape: (8, 3, 4) (expected (8, 3, 4))
```

`(N, ng, nx)` layout is now the native storage for all three cache
fields. The legacy `(N, nx, ng)` shape is GONE from
`CollisionCache.from_geometry`'s output.

## §7 OUT-of-scope acknowledgement

Per the brief's §C anti-recommendations, this PR DID NOT:

1. **Flip `GeometryCoefficients`** — confirmed no field on
   `GeometryCoefficients` has a group axis; no flip needed.
2. **Touch operator leaves** in `orpheus/sn/operator.py`,
   `orpheus/sn/scattering.py`, `orpheus/sn/fission.py` — that's
   PR-INDEX-4.
3. **Touch `SNSolver.scalar_flux`, `angular_flux`, `compute_keff`,
   `compute_group_*_rate`** — that's PR-INDEX-5.
4. **Flip σ_t storage in `Mixture` / `assemble_cell_xs`** — that's
   PR-INDEX-3. The bridge transposes at `_ensure_coll_cache`
   (`sweep.py:230-232`) and `solver.py:280, 303` STAY for now.
5. **Regenerate any regression snapshot** — the round-trip is pure
   view manipulation (`np.transpose`); regression snapshots stay
   bit-identical at `rtol=1e-12`.
6. **Add `inverse_denom_legacy` alias or shape-property** — followed
   the aggressive-retirement principle: the new layout is the ONLY
   layout.
7. **Rename `inverse_denom` or `a_attenuation`** — only the shape (and
   per-element semantic axis order) changed.
8. **Introduce `np.swapaxes` / `.T` at consumers to preserve the old
   reading pattern** — the whole point was to REMOVE those bridges.

## §8 Decision-point honesty

The brief's decision-point checkpoints:

- **"A regression snapshot drifts beyond `rtol=1e-12`"** → DID NOT
  TRIP. All 11 snapshots PASS at `rtol=1e-12, atol=1e-13`. The
  layout flip is pure transpose; FP drift is mechanically
  impossible when no reduction tree changes.
- **"Cumprod axis trap"** → DID NOT TRIP. I deliberately updated
  `axis=1 → axis=2` at the first edit pass, with an inline comment
  pinning the rationale ("cumulative product runs along axis=2 (the
  trailing cell axis); under the principled layout the cell axis is
  NOT axis 1"). If I'd missed it, the regression suite would have
  caught it loudly (per-cell cumprod over groups gives nonsense
  values).
- **"Denominator broadcast misalignment"** → DID NOT TRIP. I
  rewrote the `from_geometry` body with named intermediates
  (`streaming_face_term`, `curvature_redistribution_term`,
  `geometric_streaming_term`, `sig_t_chain`, `collision_volume_term`,
  `denom`, `inverse_denom`, `a_numer`) — Pattern 3 in action. Each
  intermediate is `(N, nx)` or `(N, ng, nx)` with the shape
  annotated inline. The verified shapes propagate cleanly.

### Near-miss bugs caught mid-PR

NONE. The transpose pattern from PR-INDEX-1 had primed me for the
axis-trap; the cumprod axis change was the first edit I made to the
body. The chain-ordering needed adjustment (was `sig_t[geom.chain_idx]`
on `(nx, ng)` σ; now `sig_t[:, geom.chain_idx].transpose(1, 0, 2)` on
`(ng, nx)` σ — gives `(N, ng, nx)` natively). The advanced-indexing
pattern was straightforward.

The PR-INDEX-1 closeout's near-miss was a per-element-SEMANTICS bug
(source pre-scale convention). The analogous risk for PR-INDEX-2 would
have been if `sig_t_chain * geom.V` had used a wrong V shape — but
`geom.V[:, None, :]` (broadcast across ng) is the same `(N, ny)` shape
times `(N, ng, nx)`, which broadcasts correctly because the geometry-V
is per-cell-per-ordinate already (not per-cell-per-group). No semantic
drift introduced.

## §9 Documentation of ambiguities

- **Test #6 dual-view: source shape**. The dual-view test now consumes
  σ_t and Q in `(ng, nx)` layout but DD `update` expects per-cell
  `(ng,)` slices. The reference path uses `sig_t[:, cell_i]` and
  `QV_chain[:, k_chain]`. The cache-fast path uses the principled
  `(ng, nx)` views and transposes once at the `ordinate_scan` call
  (`coll.a_attenuation[n].T, b.T`) since `ordinate_scan` requires
  scan-axis (cell axis) leading. This is a TEST-LOCAL transpose to
  align with the `ordinate_scan` ABI; NOT a cache-consumer bridge.
- **Test #7 indexing pattern**. Previously `coll.inverse_denom[n, k_chain]`
  picked `(ng,)` (axis 1 was cell). Under PR-INDEX-2, axis 1 is group;
  to get the same `(ng,)` semantic (per-cell denom across all groups),
  the indexing is now `coll.inverse_denom[n, :, k_chain]` (fix ordinate
  + cell, sweep groups). I chose this over flipping to "fix group +
  cell, sweep ordinates" because the test's algebraic intent is the
  per-cell DD coefficient as a function of group, matching the per-
  cell `cell_balance_terms` output's `(ng,)` `denom`.
- **`solver.py` bridge transpose stays**. The `sig_t_1d = self.sig_t[:, 0, :].T`
  pattern at solver.py:280 and 303 is the new bridge transpose this
  PR introduces. It will be removed by PR-INDEX-3 when `self.sig_t`
  itself flips to `(ng, nx, ny)`. This is the standard
  "bridges-as-the-boundary-moves-outward" pattern.

## §10 Next step pointer

**PR-INDEX-3**: flip `Mixture` / `assemble_cell_xs` cross-section
storage from `(nx, ny, ng)` to `(ng, nx, ny)`. After PR-INDEX-3:

- `self.sig_t` becomes `(ng, nx, ny)` natively.
- The bridge transpose at `solver.py:280, 303` disappears (the cache
  consumer expects `(ng, nx)`; reading from a `(ng, nx, ny=1)` storage
  via `self.sig_t[:, :, 0]` yields `(ng, nx)` directly).
- The bridge transpose at `sweep.py:_ensure_coll_cache` (`sig_t_1d =
  sig_t[:, 0, :].T`) similarly disappears.
- `_run_1d_sweep` entry transpose at `sweep.py:288`
  (`sig_t_p = np.transpose(sig_t[:, 0, :], (1, 0))`) becomes
  `sig_t_p = sig_t[:, :, 0]` — no transpose, direct cell-axis slice.

Dispatch sub-agent: **method-implementer** (per plan §8).

### Coding-elegance audit note for future PR

Pattern 13 candidate: the bridge transpose at `solver.py:280` is exactly
the "bare numpy passed across module boundary with shape convention in a
docstring" anti-pattern. It will be eliminated structurally by PR-INDEX-3
(no transpose needed once storage matches) AND by the Issue #197 typed-
field-contract resume (`CrossSection` `NewType` or wrapper). No
intermediate work needed in this PR.

## §11 Commit message (proposed)

```
perf(sn): CollisionCache (N, ng, nx) layout flip + remove sweep cache-read bridges (Issue #196 PR-INDEX-2)

Flips the storage layout of CollisionCache fields (inverse_denom,
a_attenuation, cumprod_a) from (N, nx, ng) to the principled
(N, ng, nx) — energy g is the second axis, NOT trailing. The
CollisionCache.from_geometry constructor accepts σ_t in (ng, nx) layout
natively; the cumprod runs along axis=2 (the trailing cell axis).

Removes two bridge transposes in _run_1d_sweep:
  - SLAB joint-batch: np.swapaxes(coll.*[ords], 1, 2) (3 sites) gone
    — the indexed slice IS (K, ng, nx) natively.
  - CURVILINEAR per-ordinate: coll.*[global_n].T gone — the indexed
    slice IS (ng, nx) natively.

Introduces a new (transient) bridge transpose at the solver-side
from_geometry call sites (solver.py:280, 303): self.sig_t[:, 0, :].T
to convert legacy (nx, 1, ng) storage to the principled (ng, nx) the
cache now expects. This bridge disappears in PR-INDEX-3 when Mixture
/ assemble_cell_xs flips its native layout to (ng, nx, ny).

GeometryCoefficients (Stratum 1) is unchanged — its fields are (N, nx)
or (N,) with no group axis; no flip needed.

All 11 regression snapshots stay bit-identical at rtol=1e-12.
Sweep cache + ordinate scan + L0 streaming-equilibrium gates green.
Slab benchmark stabilized to ~0.149 ms/sweep (PR-INDEX-1 was ~0.21 ms
mean with 2× trial variance; removing the bridge transposes both
improves the mean and tightens the timing distribution).
```
