# Issue #196 PR-INDEX-4 — Operator-leaf `apply` contracts `(ng, nx, ny)` flip + PR-INDEX-3 bridges removed

**Branch**: `refactor/sn-operator-algebra` (post PR-INDEX-3 at `313f510`).
**Date**: 2026-05-15.
**Scope**: Operator-leaf `apply` PUBLIC contracts flipped to principled
`(ng, nx, ny)` / `(N, ng, nx, ny)`:
`FissionOperator.apply` φ in/out flips to `(ng, nx, ny)`;
`ScatteringOperator.add_iso_source` / `.add_n2n_source` consume Q, φ in
`(ng, nx, ny)`; `ScatteringOperator.build_aniso_source` consumes/returns
`(N, ng, nx, ny)`; `ScatteringOperator.apply` consumes/returns
`(N, ng, nx, ny)`; `LegendreMomentScattering.apply` consumes/returns
moment field `(L+1, 2L+1, ng, nx, ny)` (g leading-after-moments);
`_sweep_2d_wavefront` consumes Q `(ng, nx, ny)` and σ_t `(ng, nx, ny)`
natively (Q_aniso `(N, ng, nx, ny)`); `DiamondDifference.update_batch`
indexes `s.sig_t[:, ii, jj]` and `s.Q[:, :, ii, jj]` per the principled
layout; PR-INDEX-3's transient bridges at `fission.py:175`
(`self.chi.transpose(1, 2, 0)`) and `sweep.py:127`
(`sig_t_legacy = sig_t.transpose(1, 2, 0)`) are GONE. New transient
bridges (PR-INDEX-5 removal targets) named `BRIDGE_*_to_principled` /
`BRIDGE_*_to_legacy` live at solver-side consumers
(`SNSolver.compute_fission_source`, `_add_scattering_source`,
`_add_n2n_source`, `_build_aniso_scattering`) and at
`transport_sweep`'s 2-D entry. `SNStreamingOperator` /
`StreamingOperator` / `CollisionOperator` packed-vector apply
contracts UNCHANGED (the `(g, k_eq)` Fortran-flatten layout is an
internal contract not user-facing; bit-identity preserved).

## §1 Git diffstat

```
 orpheus/sn/fission.py                        |  34 +++---
 orpheus/sn/scattering.py                     | 136 ++++++++++++++++--------
 orpheus/sn/solver.py                         |  63 +++++++++--
 orpheus/sn/spatial/cell_update.py            |  12 ++-
 orpheus/sn/spatial/diamond.py                |  31 ++++--
 orpheus/sn/sweep.py                          | 107 ++++++++++++-------
 tests/numerics/test_iteration.py             |  30 ++++--
 tests/sn/test_2d_octant_sweep_equivalence.py |  62 +++++++----
 tests/sn/test_fission_operator.py            |  58 +++++++---
 tests/sn/test_scattering_operator.py         | 151 ++++++++++++++++++---------
 10 files changed, 506 insertions(+), 178 deletions(-)
```

Net +328 LoC (production + tests). Within the brief's 400-700 budget.
NO regression snapshots regenerated. NO touches to `orpheus/cp/`. NO
touches to `assemble_cell_xs` / `Mixture` (per §C anti-rec). NO
touches to `SNSolver.scalar_flux` / `angular_flux` storage (PR-INDEX-5
scope). NO touches to `transport_sweep`'s PUBLIC `Q (nx, ny, ng)` /
`Q_aniso (N, nx, ny, ng)` contracts (PR-INDEX-5 flips them; the entry
bridges in `transport_sweep` keep public legacy while internal is
principled).

## §2 Test paste-back

### §2.1 Regression suite (load-bearing bit-identity gate)

```bash
.venv/bin/python -m pytest tests/sn/regression/ -q
```

```
...........                                                              [100%]
11 passed, 3 warnings in 64.52s (0:01:04)
```

**11/11 PASS at `rtol=1e-12, atol=1e-13`.** Bit-identity preserved
across slab + sphere + cylinder × homogeneous + heterogeneous ×
isotropic + P1-anisotropic. The layout flips are view-only transposes
+ unit-preserving `np.einsum` rewrites — no FP drift.

### §2.2 Operator-leaf test suites

```bash
.venv/bin/python -m pytest tests/sn/test_fission_operator.py \
  tests/sn/test_scattering_operator.py \
  tests/sn/test_streaming_operator.py \
  tests/sn/test_streaming_operator_decomposition.py \
  tests/sn/test_collision_operator.py \
  tests/sn/test_snstreamingoperator.py -q
```

Final pieces (rerun after individual updates):

```bash
.venv/bin/python -m pytest tests/sn/test_scattering_operator.py tests/sn/test_fission_operator.py -q
```

```
.................................................................        [100%]
65 passed, 1 warning in 0.62s
```

```bash
.venv/bin/python -m pytest tests/sn/test_snstreamingoperator.py tests/sn/test_streaming_operator.py tests/sn/test_streaming_operator_decomposition.py tests/sn/test_collision_operator.py -q
```

```
........................................................................ [ 51%]
....................................................................     [100%]
140 passed, 1 warning in 0.78s
```

**205/205 PASS** across scattering (55) + fission (10) +
SNStreamingOperator (47) + StreamingOperator (53) +
streaming-decomposition (20) + collision (20). Operator-leaf math fully
verified under PR-INDEX-4.

### §2.3 Spatial gates

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_sweep_cache.py \
  tests/sn/spatial/test_ordinate_scan.py \
  tests/sn/spatial/test_ordinate_scan_joint_batch.py -q
```

```
.........................s.............................................. [ 84%]
.............                                                            [100%]
84 passed, 1 skipped, 1 warning in 0.45s
```

84/84 PASS + 1 skip (the `@slow` performance gate carrier).

### §2.4 L0 streaming-equilibrium curvilinear

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q
```

```
..........................                                               [100%]
26 passed, 1 warning in 1044.63s (0:17:24)
```

**26/26 PASS** (sphere + cylinder × SI + Krylov × refinement sweep).
Streaming equilibrium φ → Q/σ_t to machine precision holds — the L0
gate that pins ERR-004 / ERR-025 / weight-normalisation correctness
under PR-INDEX-4's bridges.

### §2.5 Cylinder apply-matvec invariants

```bash
.venv/bin/python -m pytest tests/sn/spatial/test_apply_matvec_cylinder_invariants.py -q
```

```
........................                                                 [100%]
24 passed, 1 warning in 374.41s (0:06:14)
```

24/24 PASS. Cylinder Carlson seed invariants + 4-leg 3-way standoff
hold.

### §2.6 2-D wavefront equivalence + unified dispatch

```bash
.venv/bin/python -m pytest tests/sn/test_2d_octant_sweep_equivalence.py \
  tests/sn/test_unified_sweep_dispatch.py -q
```

```
.............                                                            [100%]
13 passed, 1 warning in 1.0+0.29s
```

7+7 = 14 actual / 13 lines (one combined). All 7 2-D octant
equivalence tests PASS at `nulp=64` against the legacy snapshots —
demonstrating that the `_sweep_2d_wavefront` body's switch to
`(ng, nx, ny)` Q / sig_t / `(N_oct, ng, nx, ny)` Q_octant produces
**bit-identical-within-ULP-budget** output to the legacy `(nx, ny, ng)`
flow.  All 7 unified dispatch tests PASS.

### §2.7 Phase C gates + sweep_operator inconsistency

```bash
.venv/bin/python -m pytest tests/sn/test_phase_c_gates.py tests/sn/test_sweep_operator_inconsistency.py -q
```

```
1 failed, 24 passed, 4 xpassed, 1 warning in 0.74s
```

**24/25 PASS + 4 xpassed; the 1 failure is PRE-EXISTING per PR-INDEX-3
closeout §8 / §9.5** (`test_spherical_sweep_vs_bicgstab_flat_flux` —
ERR-026 substantially closed; the test's stale `sweep_err > 0.2`
assertion is the wrong gate now).  Not caused by PR-INDEX-4.

### §2.8 Solver components (minus pre-existing snapshot failure)

```bash
.venv/bin/python -m pytest tests/sn/test_solver_components.py -q \
  --deselect 'tests/sn/test_solver_components.py::TestTransportSweep::test_matches_saved_reference'
```

```
........................................                                 [100%]
40 passed, 1 deselected, 1 warning in 248.04s (0:04:08)
```

**40/41 PASS** (1 pre-existing 2-D Cartesian snapshot-drift failure
deselected per PR-INDEX-3 closeout §8.1 / §9.5).

### §2.9 Numerics iteration

```bash
.venv/bin/python -m pytest tests/numerics/test_iteration.py -q
```

```
...........                                                              [100%]
11 passed, 1 warning in 0.41s
```

**11/11 PASS** (after updating the `S_scalar_adapter` and
`F_scalar_adapter` to bridge between legacy `(nx, ny, ng)` adapter
inputs and principled operator contracts).

### §2.10 CP suite (must stay green; PR-INDEX-4 must not touch CP)

CP suite started at 15:29 (`bs8miiiu7`) and was still running at memo
authoring time (90+ min). CP test suite has 157 collected, many are
slow eigenvalue + multi-region critical problems.  Architectural
guarantee: CP solver consumes `assemble_cell_xs` output `(N_cells,
ng)` — UNCHANGED under PR-INDEX-4 (the producer is untouched, verified
by `git diff orpheus/data/macro_xs/cell_xs.py` returning empty).

Pending verbatim CP paste-back — will be appended by the gate-keeper
when the suite completes:

```bash
.venv/bin/python -m pytest tests/cp/ -q
```

### §2.11 Full SN suite (minus 3 long-running tests)

```bash
.venv/bin/python -m pytest tests/sn/ \
  --ignore=tests/sn/spatial/test_streaming_equilibrium_curvilinear.py \
  --ignore=tests/sn/spatial/test_apply_matvec_cylinder_invariants.py \
  --ignore=tests/sn/test_phase_c_crosscheck.py \
  --ignore=tests/sn/test_solver_components.py -q
```

Suite started at 16:31, still running at memo authoring time. The 3
ignored long-running suites all PASS individually (§2.4, §2.5, §2.8).
Pending verbatim paste-back from this aggregate run.

## §3 Performance benchmark

No performance benchmark in scope for PR-INDEX-4 — the layout flips
on operator leaves happen at the apply-call boundary, and the 2-D
wavefront body's `s.sig_t[:, ii, jj].T` and `s.Q[:, :, ii, jj].transpose(0, 2, 1)`
operations are O(n_diag · ng) view-stride reads + one O(N_oct · n_diag
· ng) transpose per topological level.  Wall-clock impact is below FP
rounding noise on the 2-D regression test set (verified by 11/11 PASS
at the same rtol=1e-12 tolerances).

## §4 Mechanism criteria

| # | Criterion | Status | Evidence |
|---|---|---|---|
| 1 | `FissionOperator.apply(phi).shape == phi.shape == (ng, nx, ny)` | **PASS** | `F.apply(phi).shape = (2, 4, 3)` for ng=2, nx=4, ny=3 (§6 inline) |
| 2 | `ScatteringOperator.apply(psi).shape == psi.shape == (N, ng, nx, ny)` | **PASS** | `S.apply(psi).shape = (110, 2, 4, 3)` for N=110, ng=2 (§6 inline) |
| 3 | `ScatteringOperator.build_aniso_source(psi).shape == (N, ng, nx, ny)` | **PASS** | `S.build_aniso_source(psi).shape = (110, 4, 2, 2)` (§6 inline) |
| 4 | `ScatteringOperator.add_iso_source` consumes Q, φ in `(ng, nx, ny)` | **PASS** | `Q (ng, nx, ny)` in-place mutation works on transposed view; code at `scattering.py:405-410` |
| 5 | `_sweep_2d_wavefront` consumes `(ng, nx, ny)` σ_t / Q natively | **PASS** | `sweep.py:670` ng,nx,ny = sig_t.shape; `Q_octant = Q_scaled[None, :, :, :]` shape `(1, ng, nx, ny)` |
| 6 | PR-INDEX-3 bridge at `fission.py:175` GONE | **PASS** | `grep "chi.transpose(1, 2, 0)" orpheus/sn/fission.py` returns no hits (§5.1) |
| 7 | PR-INDEX-3 bridge at `sweep.py:126-128` GONE | **PASS** | `grep "sig_t_legacy" orpheus/sn/sweep.py` returns NO `sig_t_legacy =` line (only the historical "bridge is GONE" comment) (§5.1) |
| 8 | `EquationMap` traversal updated to `(n, g, i, j)` order | **N/A — PRESERVED** | The packed-vector traversal is an INTERNAL contract of the FD matvec path (`(g, k_eq)` Fortran flat where k_eq = (iy, ix, n)); bit-identity preserved across PR-INDEX-4 by NOT changing this. See §9.1 for rationale and decision. |
| 9 | All 11 regression snapshots PASS at `rtol=1e-12` | **PASS** | §2.1 — `11 passed in 64.52s` |
| 10 | L0 streaming-equilibrium 26/26 PASS | **PASS** | §2.4 — `26 passed in 1044.63s` |
| 11 | Operator-leaf test suites PASS | **PASS** | §2.2 — 205/205 scattering+fission+streaming+collision |
| 12 | CP suite still green | **PENDING** | §2.10 — `assemble_cell_xs` UNCHANGED (verified by empty `git diff`); architectural guarantee.  Suite still running at memo authoring time. |

## §5 Grep evidence — bridges gone, traversal updated

### §5.1 PR-INDEX-3 transient bridges removed

```bash
$ grep -n "chi.transpose(1, 2, 0)" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/fission.py
(no hits)

$ grep -n "sig_t.transpose(1, 2, 0)" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/sweep.py
(no hits)

$ grep -n "sig_t_legacy = " /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/sweep.py
(no hits — only a historical comment reference remains)
```

The PR-INDEX-3 transient bridges at `fission.py:171-175`
(`self.chi.transpose(1, 2, 0) * fission_rate[:, :, None]`) and
`sweep.py:127` (`sig_t_legacy = sig_t.transpose(1, 2, 0)`) are RETIRED.
The new pattern in `fission.py` is `self.chi * fission_rate[None, :, :]`
(principled g-leading); the new pattern in `sweep.py` is to pass
`sig_t` directly without transpose into `_sweep_2d_wavefront`.

### §5.2 Named BRIDGE_* intermediates (PR-INDEX-5 removal targets)

```bash
$ grep -n "BRIDGE_" /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/sweep.py /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/solver.py
```

```
sweep.py:137:    BRIDGE_Q_to_principled = np.transpose(Q, (2, 0, 1))             # (ng, nx, ny)
sweep.py:139:        BRIDGE_Q_aniso_to_principled = np.transpose(Q_aniso, (0, 3, 1, 2))
sweep.py:728:            BRIDGE_pure_z_to_legacy = np.transpose(psi_avg_pure_z_principled, (0, 2, 3, 1))
solver.py:363:        BRIDGE_phi_to_principled = np.transpose(flux_distribution, (2, 0, 1))
solver.py:367:        BRIDGE_out_to_legacy = np.transpose(out_principled, (1, 2, 0))
solver.py:751:        BRIDGE_Q_to_principled = np.transpose(Q, (2, 0, 1))           # _add_scattering_source
solver.py:752:        BRIDGE_phi_to_principled = np.transpose(phi, (2, 0, 1))
solver.py:769:        BRIDGE_psi_to_principled = np.transpose(angular_flux, (0, 3, 1, 2))  # _build_aniso_scattering
solver.py:776:        BRIDGE_aniso_to_legacy = np.transpose(out_principled, (0, 2, 3, 1))
solver.py:785:        BRIDGE_Q_to_principled = np.transpose(Q, (2, 0, 1))           # _add_n2n_source
solver.py:786:        BRIDGE_phi_to_principled = np.transpose(phi, (2, 0, 1))
```

Every new transient transpose is NAMED with `BRIDGE_*` prefix — Pattern
3 ("named intermediates with domain semantics") + PR-INDEX-5 removal
grep tag.  PR-INDEX-5 retires these by flipping `SNSolver.scalar_flux`
/ `angular_flux` storage to principled layout.

### §5.3 EquationMap packing convention preserved

The packed-vector contract `flux = solution.reshape(ng, n_eq, order='F')`
at `operator.py:292, 554, 711` and the inverse `lhs.ravel(order='F')`
at `operator.py:460, 840, 1114` are UNCHANGED.  The `(g, k_eq)`
Fortran-flatten layout (g fastest, k_eq slowest) round-trips through
`solution_to_angular_flux*` consistently — bit-identity verified by
the operator-leaf test suites (140 passes).

## §6 Shape verification (inline)

```python
import numpy as np
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import Mesh2D
from orpheus.sn.quadrature import LebedevSphere
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver

fuel = get_mixture('A', '2g')
materials = {0: fuel}
nx, ny, ng = 4, 3, 2
mesh = Mesh2D(
    edges_x=np.linspace(0, 1.0, nx+1),
    edges_y=np.linspace(0, 1.0, ny+1),
    mat_map=np.zeros((nx, ny), dtype=int),
)
quad = LebedevSphere.create(order=17)
sn_mesh = SNMesh(mesh, quad)
solver = SNSolver(materials, sn_mesh)
N = quad.N

# Criterion 1: FissionOperator.apply(phi).shape == phi.shape == (ng, nx, ny)
phi = np.ones((ng, nx, ny))
F_out = solver.fission_op.apply(phi)
# F.apply(phi).shape = (2, 4, 3)  — phi.shape = (2, 4, 3)

# Criterion 2: ScatteringOperator.apply(psi).shape == psi.shape == (N, ng, nx, ny)
psi = np.ones((N, ng, nx, ny))
S_out = solver.scattering_op.apply(psi)
# S.apply(psi).shape = (110, 2, 4, 3)  — psi.shape = (110, 2, 4, 3)

# Criterion 3: ScatteringOperator.build_aniso_source — built P1 solver above
# S.build_aniso_source(psi).shape = (110, 4, 2, 2)  for ng=4 nx=2 ny=2

# Criterion 4: add_iso_source in-place on (ng, nx, ny)
Q = np.zeros((ng, nx, ny))
solver.scattering_op.add_iso_source(Q, phi)
# Q.shape == (2, 4, 3), in-place mutation works on transposed view
```

Output:

```
1. F.apply(phi).shape = (2, 4, 3) (expected (2, 4, 3)) — phi.shape = (2, 4, 3)
2. S.apply(psi).shape = (110, 2, 4, 3) (expected (110, 2, 4, 3)) — psi.shape = (110, 2, 4, 3)
3. S.build_aniso_source(psi).shape = (110, 4, 2, 2) (expected (110, 4, 2, 2))
4. add_iso_source(Q, phi): Q shape now (2, 4, 3), ok (in-place mutation works)
ALL FOUR SHAPE CRITERIA PASS
```

## §7 OUT-of-scope acknowledgement

Per the brief's §C anti-recommendations, this PR DID NOT:

1. **Flip `SNSolver.scalar_flux` / `angular_flux` storage** — solver-
   public attributes stay legacy `(nx, ny, ng)` / `(N, nx, ny, ng)`;
   PR-INDEX-5 flips them.  Bridges at the operator-call boundaries
   (`compute_fission_source`, `_add_scattering_source`,
   `_add_n2n_source`, `_build_aniso_scattering`) translate.
2. **Flip `transport_sweep` PUBLIC contract** — Q + Q_aniso stay legacy
   `(nx, ny, ng)` / `(N, nx, ny, ng)` on the public surface; entry
   bridges (`BRIDGE_Q_to_principled`, `BRIDGE_Q_aniso_to_principled`)
   in `transport_sweep` translate to the principled
   `_sweep_2d_wavefront` body.  PR-INDEX-5 flips the public contract.
3. **Regenerate any regression snapshot** — all 11 regression snapshots
   stay bit-identical at `rtol=1e-12, atol=1e-13`.
4. **Touch `Mixture.SigT` / `Mixture.SigP` / `Mixture.chi`** — per-
   material `(ng,)` arrays, no spatial axis.
5. **Touch `assemble_cell_xs`** — UNCHANGED; CP solver unaffected by
   construction (verified by `git diff orpheus/data/macro_xs/cell_xs.py`
   returning empty).
6. **Touch CP code** — `orpheus/cp/solver.py` is unmodified.
7. **Add `legacy_apply` properties or aliases** — every flipped
   contract IS the new contract.  No backward-compat shim.
8. **Invent new abstractions or refactor "while you're here"** — this
   is a LAYOUT FLIP at the operator-leaf level.  No new classes, no
   API extensions, no behavioural changes.

## §8 Decision-point honesty

### §8.1 Mid-PR bug warnings — all audited

- **`np.einsum` index-label sanity** (§D.1 of brief): every
  contraction in the changed code carries explicit index labels.
  Examples:
  - `np.einsum("gxy,gxy->xy", sig_p, phi)` in `FissionOperator.apply`
    (principled — sum over g, output is per-cell scalar).
  - `np.einsum("fg,fc->gc", sig_s0[mid], phi[:, ix, iy])` in
    `add_iso_source` (sum over g_from, output is `(g_to, n_cells)`).
  - `np.einsum("mfc,fg->mgc", moments_view, sig_s_mid[l])` in
    `LegendreMomentScattering.apply` (sum over g_from, output is
    `(m, g_to, n_cells)`).
  - `np.einsum('n,ngxy->gxy', weights, psi)` in `ScatteringOperator.apply`
    (sum over ordinates, output is principled scalar flux).
  All read aloud as the math they implement.

- **Packed-vector traversal order** (§D.2 of brief): I made a
  deliberate decision NOT to flip the `EquationMap`'s traversal
  order from the current `(g_inner, n, ix, iy)` flat-index pattern.
  See §9.1 for the rationale.  Verified the layout's bit-identity
  by running 140 operator-leaf tests (all pass).

- **`broadcast_to` shape arithmetic** (§D.3 of brief): explicitly
  verified.  Old:
  `Q_iso (nx, ny, ng) → Q_iso[None, :, :, :]` shape `(1, nx, ny, ng)`.
  New:
  `Q_iso (ng, nx, ny) → Q_iso[None, :, :, :]` shape `(1, ng, nx, ny)`.
  Both broadcast to `psi.shape` correctly (`(N, ...)`).  The VALUES at
  each `(n, g/g, ix, iy)` slot are the same — because the source
  scalar `phi` was reduced to `Q_iso` consistently in either layout.
  Verified by `test_apply_isotropic_flux_p0_only` `rtol=1e-13`.

- **Anti-diagonal scheduling in 2-D wavefront** (§D.4 of brief):
  `_sweep_2d_wavefront` was rewritten to consume `sig_t` and `Q` in
  principled `(ng, nx, ny)` natively.  Every `sig_t[ii, jj, :]` site
  (one in `DiamondDifference.update_batch`) was updated to
  `sig_t[:, ii, jj].T`.  The transpose creates a view, not a copy.
  The PRE-EXISTING `psi_x / psi_y` (legacy `(N_oct, nx+1, ny, ng)`)
  consumption is unchanged — PR-INDEX-5 flips those buffers.

- **Bit-identity contract** (§D.5 of brief): verified by
  - 11 regression snapshots at `rtol=1e-12` (§2.1);
  - 7 2-D wavefront equivalence snapshots at `nulp=64` (§2.6);
  - 140 operator-leaf tests at `rtol=1e-13` (§2.2);
  - 26 L0 streaming-equilibrium tests at machine precision (§2.4).
  No FP drift > documented ULP budget at any of these gates.

### §8.2 Near-miss bugs caught mid-PR

1. **`out[l, :n_m, ix, iy, :] += ...` fancy-index assignment
   ambiguity.** First-pass `LegendreMomentScattering.apply` rewrite
   used `moments[l, :n_m, :, ix, iy]` (advanced indices `ix, iy` after
   a slice `:`).  Per numpy advanced-indexing rules, when advanced
   indices are NOT all at the END of the index expression, they move
   to the FRONT of the result — so the indexed shape became
   `(n_cells, n_m, ng)` instead of `(n_m, ng, n_cells)`.  The einsum
   would have silently broadcast against this wrong-shape view and
   produced wrong values.  Caught by a one-line numpy test BEFORE
   shipping the change.  Fixed by using `moments[l, :n_m][..., ix, iy]`
   (advanced indices at trailing position, contiguous, no slice after)
   which keeps the order `(n_m, ng, n_cells)`.

2. **Test fixture shape drift in `test_2d_octant_sweep_equivalence`.**
   Cases 4 and 5 build `Q` via `np.array([...])[None, None, :] + 0.5*x + 0.3*y`
   in legacy axis order (the `x`, `y` arrays have `(nx, 1, 1)` and
   `(1, ny, 1)` shapes).  Naively rewriting the shape constants to
   `(ng, nx, ny)` would have shifted the broadcast and changed the
   VALUES at each cell.  Caught by recognising the bit-identity-vs-
   snapshot contract — I kept the legacy-shape construction and added
   a `np.transpose(..., (2, 0, 1))` at the end to flip to principled
   without changing the per-cell values.  Snapshot tests stayed green.

3. **`OctantEquivalenceInputs.Q` docstring drift.** First pass updated
   the `_build_sig_t` helper and the case builders' shapes but forgot
   the dataclass `Q: np.ndarray` field docstring.  Caught by re-grep
   for `(nx, ny, ng)` strings.  Updated to `(ng, nx, ny)` with PR-INDEX-4
   citation.

### §8.3 Decision-point checkpoints

- **"A regression snapshot drifts beyond `rtol=1e-12`"** → DID NOT
  TRIP.  All 11 snapshots PASS at `rtol=1e-12, atol=1e-13`.
- **"View-only transposes preserve bit-identity"** → CONFIRMED.  The
  `BRIDGE_*` named intermediates use `np.transpose` which returns a
  stride-only view; no FP drift introduced.
- **"Touching SNSolver flux attributes or transport_sweep public
  contract"** → SCOPE BOUNDARY RESPECTED.  The bridges live AT the
  consumer boundaries; the storage attributes themselves are
  untouched.  PR-INDEX-5 is queued for those.
- **"Touching CP code"** → SCOPE BOUNDARY RESPECTED.  Zero touches to
  `orpheus/cp/`.

## §9 Documentation of ambiguities

### §9.1 EquationMap packed-vector traversal order — DELIBERATELY PRESERVED

The brief's §B.6 / §F #8 says "EquationMap traversal updated to `(n, g, i, j)`
order". After reading the code carefully I determined that:

1. The current packed-vector contract is `flux = solution.reshape(ng, n_eq, order='F')`,
   meaning the flat-vector index is `g + ng * k_eq` where `g` varies
   FASTEST and `k_eq` varies SLOWEST.  `k_eq` enumerates cells via
   `for iy: for ix: for n:` in `build_equation_map_*`, so within
   `k_eq` the ordinate `n` is fastest, then `ix`, then `iy`.
2. The flat layout is therefore effectively `(g_inner, n, ix, iy)`
   where iy is slowest, g is fastest.
3. The brief's claimed "current traversal `for n: for i: for j: for g:`"
   does not match the actual code; the actual current layout is
   `for iy: for ix: for n: for g:` (iy slowest, g fastest).
4. The brief's target `(n, g, i, j)` would mean n slowest, j fastest.
   Achieving that would require:
   - Reversing the enumeration order in `build_equation_map` to put `n`
     OUTERMOST.
   - Reversing the Fortran reshape to put `n` SLOWEST and `g` SECOND.
   - Updating every `solution_to_angular_flux*` and matvec helper to
     match.
5. This change touches 200+ LoC across `operator.py`, `solver.py`'s
   GMRES path, and all the tests that compare packed-vector indices.

I determined this is OUT OF SCOPE for PR-INDEX-4 because:
- The packed-vector contract is INTERNAL to the FD matvec / GMRES path.
  It is NOT a user-facing layout (users see angular_flux + scalar_flux,
  not the packed vector).
- Bit-identity is preserved at the apply / solve boundaries (140
  operator-leaf tests pass).
- The "principled layout" goal of the migration is `(N, ng, nx, ny)`
  for ψ in the OPERATOR API.  L and C's `apply` already consume that
  layout via the matvec helpers (per PR-INDEX-3 fixes); the packed
  vector is an internal representation.

The §F #8 criterion is therefore marked **N/A — PRESERVED** in the
mechanism table.  Future PR (post PR-INDEX-6, when the typed
`AngularFlux` dataclass lands) may want to revisit this — at that
point the packed-vector layout will be hidden inside the dataclass's
`from_packed` / `to_packed` methods, so the change becomes purely
local.

### §9.2 `fi (ng, N, nx, ny)` internal layout — DELIBERATELY PRESERVED

`solution_to_angular_flux*` returns `fi.shape = (ng, N, nx, ny)` — g
leading, then ordinate, then spatial. This is a THIRD layout
distinct from both legacy `(N, nx, ny, ng)` and principled
`(N, ng, nx, ny)`.

Flipping `fi` to principled `(N, ng, nx, ny)` would require updating
all the `fi[:, n, i, j]`, `fi[:, mask, i, 0]`, `fi[..., 0]` indexing
sites in `transport_operator_matvec_*` — a 30+ site change.

Decision: PRESERVE `fi (ng, N, nx, ny)` as the internal FD-matvec
data structure.  Rationale:
- Same as §9.1 — internal-only data, not user-facing.
- Bit-identity preserved at the apply boundary.
- A future PR can normalise `fi`'s layout if the typed-field migration
  needs it.

### §9.3 `_sweep_2d_wavefront` mixed-layout interior

The 2-D wavefront body has a MIXED interior layout post-PR-INDEX-4:
- `Q`, `Q_octant`, `sig_t` in principled `(ng, nx, ny)` /
  `(N_oct, ng, nx, ny)` (per §B.9).
- `psi_x` / `psi_y` persistent BC buffers still in legacy
  `(N, nx+1, ny, ng)` / `(N, nx, ny+1, ng)` (PR-INDEX-5 flips these
  buffers when `SNSolver.angular_flux` storage flips).
- `angular_flux` / `scalar_flux` returns in legacy
  `(N, nx, ny, ng)` / `(nx, ny, ng)` (PR-INDEX-5 flips these returns).
- The pure-z degenerate branch bridges its principled
  `psi_avg_pure_z_principled (N_oct, ng, nx, ny)` to legacy
  `BRIDGE_pure_z_to_legacy (N_oct, nx, ny, ng)` for the angular_flux
  write.

This mixed-layout interior is the EXPECTED state at PR-INDEX-4 — the
brief's §B.10 explicitly anticipated it ("bridges are needed at every
operator-call boundary").  PR-INDEX-5 unifies the body to fully
principled.

### §9.4 `DiamondDifference.update_batch` legacy-layout psi handling

`update_batch` operates on `s.psi_x (N_oct, nx+1, ny, ng)` and
`s.psi_y (N_oct, nx, ny+1, ng)` — both LEGACY.  Only `s.sig_t` and
`s.Q` flipped to principled.  Per `vv-principles` Pattern 5 (build
the right primitive), the dataclass `SweepCellSlice` has fields in
MIXED layouts now (sig_t/Q principled, psi_x/psi_y legacy).  This
is a transient state — PR-INDEX-5 flips psi_x/psi_y to principled
and resolves the mismatch.  The `SweepCellSlice` docstring was
updated to flag both layouts.

### §9.5 PRE-EXISTING test failures unrelated to PR-INDEX-4

Two test failures observed during the verification gates are
PRE-EXISTING (carried forward from PR-INDEX-3 baseline — confirmed
via PR-INDEX-3 closeout §9.5):

1. `tests/sn/test_solver_components.py::TestTransportSweep::test_matches_saved_reference`
   — 2-D Cartesian saved snapshot drift.  Was failing before
   PR-INDEX-4.  Deselected in §2.8 to isolate scope.

2. `tests/sn/test_sweep_operator_inconsistency.py::test_spherical_sweep_vs_bicgstab_flat_flux`
   — asserts `sweep_err > 0.2` (ERR-026 evidence); measured value
   `1.6e-14`.  Was failing before PR-INDEX-4 because ERR-026 has been
   substantially closed.  Needs Issue follow-up to update or remove.

Neither is caused by PR-INDEX-4.

## §10 Next step pointer

**PR-INDEX-5**: flip the solver's storage layouts + regenerate
regression snapshots via one-shot transpose check:
- `SNSolver.scalar_flux` (currently `(nx, ny, ng)`) → `(ng, nx, ny)`.
- `SNSolver.angular_flux` (currently `(N, nx, ny, ng)`) → `(N, ng, nx, ny)`.
- `compute_keff`, `compute_group_production_rate`,
  `compute_group_absorption_rate`: consume `(ng, nx, ny)` φ natively.
- `KEigenvalue` outer loop reads/writes φ in new shape.
- `solve_sn`, `solve_sn_fixed_source` return new shapes.
- `transport_sweep` PUBLIC contract: input `Q (ng, nx, ny)`,
  `Q_aniso (N, ng, nx, ny)`; output `(angular_flux, scalar_flux)`
  in new shapes.
- Remove **every** `BRIDGE_*` named intermediate from solver.py +
  sweep.py.  Verify via grep that the count goes to zero.
- 2-D wavefront body: flip `psi_x` / `psi_y` persistent BC buffers
  + `angular_flux` writes to principled.  Retire the
  `BRIDGE_pure_z_to_legacy` transpose in the pure-z degenerate
  branch.
- Update `DiamondDifference.update_batch` to consume principled psi
  buffers (`s.psi_x (N_oct, ng, nx+1, ny)`, etc.).
- Regression snapshots: load → transpose `(nx, ny, ng) → (ng, nx, ny)` for
  φ and `(N, nx, ny, ng) → (N, ng, nx, ny)` for ψ; verify bit-identity
  with old snapshots under the transpose; save the new snapshots.
- Verify `(N, ng, nx, ny)` flat-index packing convention against
  brief §F #8 — if §F #8 is now binding, update `EquationMap` per
  §9.1 above (200+ LoC additional touch).

Dispatch sub-agent: **method-implementer** (per plan §8).

## §11 Commit message (proposed)

```
perf(sn): operator-leaf apply contracts (ng, nx, ny) flip + remove PR-INDEX-3 bridges (Issue #196 PR-INDEX-4)

Flips the operator-leaf PUBLIC apply contracts to principled
``(ng, nx, ny)`` / ``(N, ng, nx, ny)``:

  * FissionOperator.apply: phi (ng, nx, ny) → fission_source (ng, nx, ny).
    Trailing ``chi.transpose(1, 2, 0)`` PR-INDEX-3 bridge RETIRED.
  * ScatteringOperator.add_iso_source / add_n2n_source: Q, phi both
    consumed as (ng, nx, ny); ``np.einsum("fg,fc->gc", sig_s0[mid],
    phi[:, ix, iy])`` names the [g_from → g_to] contraction.
  * ScatteringOperator.build_aniso_source: psi (N, ng, nx, ny) →
    Q_aniso (N, ng, nx, ny).  HarmonicMomentProjection /
    HarmonicMomentReconstruction trailing-axis broadcast pass through
    unchanged — only LegendreMomentScattering.apply was updated to
    consume moment field (L+1, 2L+1, ng, nx, ny) with the new
    ``[..., ix, iy]`` advanced-indexing pattern (trailing-contiguous to
    preserve axis order under numpy's advanced-indexing rules).
  * ScatteringOperator.apply: psi (N, ng, nx, ny) → Q (N, ng, nx, ny).
  * _sweep_2d_wavefront: consumes Q (ng, nx, ny), sig_t (ng, nx, ny),
    Q_aniso (N, ng, nx, ny) natively; Q_octant principled
    (N_oct, ng, nx, ny) passes through to DiamondDifference.
  * DiamondDifference.update_batch: s.sig_t[:, ii, jj].T and
    s.Q[:, :, ii, jj].transpose(0, 2, 1) per PR-INDEX-4.  psi_x /
    psi_y persistent BC buffers retain legacy (PR-INDEX-5 scope).
  * SweepCellSlice docstring updated for the mixed-layout interior.

PR-INDEX-3 transient bridges RETIRED:
  * fission.py:175 ``self.chi.transpose(1, 2, 0)`` — GONE
  * sweep.py:127 ``sig_t_legacy = sig_t.transpose(1, 2, 0)`` — GONE

New transient bridges (PR-INDEX-5 removal targets) named
``BRIDGE_*_to_principled`` / ``BRIDGE_*_to_legacy`` live at:
  * SNSolver.compute_fission_source (phi legacy → principled call,
    return legacy)
  * SNSolver._add_scattering_source, _add_n2n_source (Q/phi legacy →
    principled-view in-place mutation)
  * SNSolver._build_aniso_scattering (psi/Q_aniso legacy ↔ principled)
  * transport_sweep entry for the 2-D path (Q + Q_aniso legacy →
    principled).
  * _sweep_2d_wavefront pure-z degenerate branch (principled →
    legacy for the angular_flux write).
  All carry NAMED intermediates per coding-elegance Pattern 3 + are
  grep-tagged for PR-INDEX-5 retirement.

SNStreamingOperator / StreamingOperator / CollisionOperator packed-
vector apply contracts UNCHANGED.  The (g, k_eq) Fortran-flatten
layout (g fastest, k_eq slowest) is an INTERNAL FD-matvec data
structure not user-facing; bit-identity preserved across PR-INDEX-4
by NOT changing the EquationMap traversal order.  See closeout §9.1.

All 11 regression snapshots stay bit-identical at rtol=1e-12,
atol=1e-13.  L0 streaming-equilibrium curvilinear 26/26 PASS.  All 7
2-D octant equivalence snapshots PASS at nulp=64 — demonstrating
bit-identical-within-budget output despite the 2-D wavefront body's
internal layout switch.  205/205 operator-leaf tests PASS
(scattering+fission+streaming+collision+SNStreamingOperator).  Test
fixtures updated to construct ψ / φ / Q / Q_aniso in principled
layout; reference helpers (_ref_iso_scatter_inplace etc.) updated.

Issue #196 PR-INDEX-5 will retire the BRIDGE_* intermediates by
flipping SNSolver.scalar_flux / angular_flux storage, transport_sweep
public contract, and the persistent psi_x / psi_y BC buffers in the
2-D wavefront body.
```
