# N-D layout foundation — `(N, ng, *spatial)` array convention (eliminate the phantom `ny=1`)

**Branch/worktree:** `worktree-sn-nd-layout` (fresh off `main` @ `f0ffb30`).
**Mode:** main-agent direct authorship, turn-by-turn user steering
(`feedback_no_method_implementer_for_surgical_carves`). Commit per stage; gates
throughout. Standing commit exclusions: `.claude/agent-memory/`, `.mcp.json`,
hooks, `derivations/diagnostics/`.

## Goal (user, 2026-06-09)
Make EVERY array genuinely `(N, ng, *spatial)`, spatial rank == ndim: `(nx,)`
1-D, `(nx,ny)` 2-D, `(nx,ny,nz)` future 3-D. **NOT adding 3-D compute** — a
dimension-agnostic FRAMEWORK that takes 3-D later. **Full rank-d, all paths incl.
curvilinear** ("how we get to a good architecture").

## The crux (empirically confirmed)
The pack convention (`SNMesh.spatial_shape`/`n_unknowns_flat`, from the ultraplan
C1 already on `main`) is genuine rank-d, but the typed FIELD arrays
(`_bases.py:189/255 _shape_for_mesh` → `(N,ng,nx,ny)`) are ALWAYS rank-2: 1-D
carries a phantom `ny=1` (`AngularFlux.values.shape==(4,1,5,1)` for a 5-cell
slab). **2-D is already `(N,ng,*spatial)`-correct → stays bit-identical.** The
phantom lives entirely on the 1-D path; curvilinear is also `ndim=1` → slab +
curvilinear move to rank-d TOGETHER (shared `_shape_for_mesh`).

## Operator-machinery audit (explorer, the partition the staging rests on)
**THE LEVER:** `_bases.py` `BulkField._shape_for_mesh`→`(N,ng,*spatial_shape)`,
`ScalarField._shape_for_mesh`→`(ng,*spatial_shape)`. `Field.__post_init__` enforces
`values.shape==space.shape` STRICTLY → flipping the lever makes every stale rank-2
producer **fail LOUD**. No load-bearing `ny=1` broadcast exists — the phantom is
vestigial positional indexing.

**SURVIVE UNTOUCHED (rank-agnostic already):** `FissionOperator`+`RankOneOperator`
(reduce axis=0, broadcast *spatial); all `MomentProjection`/`HarmonicMomentReconstruction`
einsums (ellipsis `n...`); `ordinate_scan`/`_pair_monoid_scan`; the
`MorelMontryAngularSweep` pole closure INTERNALS (already rank-1 — drops y on
entry: `psi_view.transpose(1,0,2,3)[...,0]`); all iteration drivers
(`SourceIteration`/`KrylovAcceleration`/`KEigenvalue`, via the ravellable
protocol); `TimedFullField.to_flat`/`from_flat` (ravel/reshape vs template);
the 1-D face layout (already `(N,ng)`); the ENTIRE 2-D Cartesian path.

**MUST CHANGE:**
- *Rank-ADAPTIVE-safe (works for current rank-2 AND future rank-1, bit-identical
  on both → prep commit C1):* `operator.py:1963/1994` CollisionOperator
  `sigma[None,:,:,:]`→`sigma[None]`; `angular_flux.py:131` `"n,ngij->gij"`→
  `"n,ng...->g..."`; any sibling `[None,:,:,:]`/`ngij`/`gxy` in scattering/sources.
- *Coupled to the lever (sigma rank must match ψ rank → C2):*
  `material_xs_field.py:369-375` `(ng,nx,ny)` reshape → `(ng,*spatial_shape)` +
  the per-material `(ix,iy)` verbs + `cells_by_material` → ndim-generic
  (`np.nonzero(mat_map)` → ndim index tuple); `geometry.py __init__`
  `(nx,ny)`/`ny=1` normalization (`volumes`/`mat_map` reshape, `is_1d`).
- *Structured 1-D (the rank change body → C2):* the `operator.py` 1-D matvec
  triple `_compute_LpC`/`_compute_decomposition`/`_compute_LpC_transpose`
  (`np.zeros((ng,N,nx,ny))` ny=1, `psi_view[:,:,0,0]` pole seed →`[:,:,0]`,
  `[:,global_dir,i,0]` writes, `sigma_t[:,:,0]`, `volumes[:,0]`,
  `.transpose(1,0,2,3)[...,0]`, re-pad `[...,None]`); `sweep.py
  _sweep_1d_unified` entry/exit (`[:,:,:,0]` drop + `[:,:,None]` re-add).
- *Solver (C2):* `solver.py` initial-guess `np.ones((ng,nx,ny))`, shape validators
  `expected=(N,ng,nx,ny)`, rate reductions `reshape(nx*ny,ng)`→`prod(spatial_shape)`.
- *Tests (C2):* slab/sphere/cylinder tests hardcoding `(N,ng,nx,1)`.

**Latent (NOT the phantom, NOT load-bearing, leave for a separate cleanup):**
`to_flat.size != n_unknowns_flat` (boundary block: `n_unknowns_flat` counts
outflow-only ordinates, `BoundaryFlux` stores all N). `n_unknowns_flat` has ZERO
production callers. Bulk parts coincide → rank change neither causes nor fixes it.

## Staging (commits, each gated)
- **C1 — local-operator rank-adaptive prep** (bit-identical on 2-D AND current
  rank-2 1-D; ALL tests green). `sigma[None]`, ellipsis einsums. Checkpoint.
- **C2 — THE LEVER + coupled + structured 1-D + tests** (the rank change). Flip
  `_shape_for_mesh`; fix material_xs, geometry `__init__`, the 1-D matvec triple,
  `_sweep_1d_unified`, solver factories/validators; migrate 1-D tests. Curvilinear
  sequenced LAST within C2, most-gated. 2-D bit-identical; 1-D re-validated via
  `k_inf=νΣ_f/Σ_a` + MMS (NOT old-vs-new bytes).
- **C3 — sweep DAG d-generic** (`OctantLabel(signs:tuple)`, `from_cartesian_2d`→
  `from_cartesian`, `diamond.cell_kernel_batch` `denom=σt+Σ_a s[a]`,
  `SweepCellSlice` per-axis). 2-D bit-identical (sum in x,y order); 3-D shape-admit.
- **C4 — boundary inventory #220** (`_resolve_bcs`→`dict[FaceLabel,op]` loop,
  `boundary_face_layout`→`face_labels`; retire named `bc_xmin/...`).
- **C5 — 3-D shape-admission pins + docs** (synthetic-3-D pack/level pins, no
  compute; extend `docs/theory/index_convention.rst` `(N,ng,nx,ny)`→`(N,ng,*spatial)`).

## Gates (test-architect, paths confirmed on `main`)
- **2-D bit-identity moat (every commit):** `tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py`
  (nulp=64) + `tests/sn/regression/test_dd_regression.py` `-W error::DriftWarning`.
- **Sentinels WITHOUT `-O`** (Mode 8): `pytest -m sentinel`. Incl. per-ordinate
  flat-flux `(L+C)·1=σ_t·1` (ERR-026), closed-form `test_kinf_homogeneous.py`.
- **1-D re-validation (C2):** `k_inf=νΣ_f/Σ_a` (`tests/sn/eigenvalue/`), MMS
  (`tests/sn/verification/mms/` slab+sphere+cyl), SI≡Krylov twin.
- **Fast per-commit (`-O`):** `tests/sn/primitives/ tests/sn/sweep/core/
  tests/sn/sweep/slab/ tests/sn/sweep/curvilinear/ tests/sn/solve/`.
- **L16 wall-clock:** the `Σ_axis` stays a `d≤3` vectorized loop, NO per-cell Python.
- Non-square `(5,7)` for Mode-6 `order=` pins; extent-1 `(1,7)`/`(nx,1)` for #214.
- Default `python -O`; env `PYTHONPATH=$worktree .venv/bin/python`. NO
  `continuous_get` (#212). Standing reds (don't block, confirm red at clean HEAD):
  #206 cyl-matvec, #195 MMS@160, #212, #214.

## References
- `.claude/plans/nd_foundation.md` (Phase-6 design-of-record §1–§2),
  `.claude/plans/sn_development_sequence.md` (Phase 6 + verification spine).
- explorer memo `sn_phantom_axis_rank_change_audit.md`; `[[lessons-L17]]`
  (crosswalk before carve), `[[lessons-L24]]` (re-characterise at pickup),
  `[[feedback-aggressive-retirement]]`.
