---
name: scan-march-seam-crosswalk
description: Structural seam crosswalk for the ScanMarch SweepStrategy (#222). The DD-kernel↔scan-recurrence α/β algebra (α=2s_x/D−1, β=2(Q+s_y·ψy_in)/D, verified at nulp); a_attenuation cache is 1-D-ONLY today (2-D recomputes inline, #206); the march+reflective-shed seam vs the wavefront capture_x/capture_y; the strategy slot; and the anti-surprises (psi_avg=0.5(in+out) STILL holds in 2-D; anisotropy/moment/multigroup are OUTSIDE the recurrence; octant=independent; curvilinear is d=1-only).
metadata:
  type: project
---

Seam crosswalk for building the **`ScanMarch`** SweepStrategy (issue #222, the S5
headline of `.claude/plans/sn_sweep_strategy.md`). Read-only investigation,
worktree `worktree-sn-nd-layout`, 2026-06-10. COMPLEMENTS the test-architect's
`scan_march_verification.md` (the V&V gate table) — this is the STRUCTURAL build
crosswalk. All `file:line` verified against the worktree; the SHAPE is durable,
line numbers drift (re-confirm via Nexus `context`).

## ⭐ EXECUTIVE SUMMARY (the load-bearing results)

1. **α/β crosswalk — CONFIRMED EXACTLY** (numerically, modulo FP-association →
   principled-equiv-not-bit-identity, the right gate). For the per-octant x-scan
   at a fixed y-row, with `D = σ_t + s_x + s_y` (the 2-D kernel denominator):
   - **α(i) = 2·s_x(i)/D(i) − 1** ✅ (= `a_attenuation` at `s_y=0`; degenerates to
     the slab cache `a = 2s/(s+σ_t)−1` exactly when the transverse axis is absent).
   - **β(i) = 2·(Q(i) + s_y(i)·ψy_in(i)) / D(i)** ✅ — the transverse-y inflow
     coupling `s_y·ψy_in` is folded into the affine SOURCE exactly the way the 1-D
     curvilinear path folds the M-M angular contribution into `b`
     (`b = 2·(QV + ang_contrib)·inv_denom`). The scan recurrence
     `out_x(i+1) = α(i)·out_x(i) + β(i)`, seeded `out_x(0)=ψ_x_in_boundary`, is
     `ordinate_scan(α, β, ψ_0)` — BYTE-for-byte the 1-D primitive.
2. **`a_attenuation` cache is 1-D-ONLY today.** The two-stratum cache
   (`CollisionCache.from_geometry`, `sweep_cache.py:438-449`) is built and consumed
   exclusively by `CumprodScan`/`_sweep_1d_unified`; the 2-D wavefront recomputes
   `D`/`a` inline per cell per level per sweep (`apply_windowed`'s
   `cell_kernel_batch` call, `sweep_graph.py:905-916`). This inline recompute IS the
   #206 duplication. The scan-march adopts the cache by construction — but it needs a
   **2-D extension** (`s_y·ψy_in` makes `D` depend on BOTH `s_x` and `s_y`, and `α`
   on both).
3. **Top-3 anti-surprises:** (a) `psi_avg = 0.5·(in_x + out_x)` STILL holds in 2-D
   (trivially, from the WDD closure) — the transverse coupling lives entirely in β,
   NOT in the averaging; (b) anisotropy/scattering/fission/multigroup/moment-output
   are ALL outside the recurrence (Q is pre-assembled, moments are an OUTPUT
   projection) — the scan-march handles the SAME per-ordinate scalar recurrence the
   wavefront does; (c) the y-out face `out_y = 2·psi_avg − ψy_in` must be threaded
   row→row in the octant's y-direction, and reflective-BC shedding moves from the
   wavefront's per-level `capture_x`/`capture_y` to "last x-scan value + last y-row
   ψy" (ERR-056 shed-order is a d≥3-only stressor in pure 2-D Cartesian).

---

## 1. DD-KERNEL → SCAN-RECURRENCE CROSSWALK (the artifact)

### 1.1 The 2-D DD cell kernel (the GROUND TRUTH)

`DiamondDifference.cell_kernel_batch` (`orpheus/sn/spatial/diamond.py:273-329`),
d-generic, but at d=2 (per axis `a∈{x,y}`):

```
denom = sigt_cells                               # (ng, n_diag)   — Σ_t on the level
numer = Q_cells                                  # (N_oct or 1, ng, n_diag)
for s_a, in_a in zip(s_axes, psi_in):            # LEFT FOLD (bit-id load-bearing)
    denom = denom + s_a                          # → D = σ_t + s_x + s_y
    numer = numer + s_a * in_a                   # → Q + s_x·in_x + s_y·in_y
psi_avg = numer / denom                          # ψ̄ = (Q + s_x·in_x + s_y·in_y)/D
psi_out = tuple(2.0*psi_avg - in_a for in_a in psi_in)   # out_a = 2ψ̄ − in_a
```

with `s_axes[a] = 2|μ_a|/Δ_a` = `streaming_x`/`streaming_y` sliced to the octant's
ordinates (`geometry.py:1347-1354`; the `s_a = 2|μ|/Δ_a` accessor is
`SNMesh.streaming(axis)` at `geometry.py:660-687`, wrapping `streaming_x`/`y`).
Terms + shapes:

| term | shape | meaning |
|------|-------|---------|
| `D = σ_t + s_x + s_y` | `(N_oct, ng, n_diag)` | denominator; the LEFT-FOLD order `((σ_t+s_x)+s_y)` is bit-id load-bearing (do NOT `sum()`) |
| `s_x = 2|μ_x|/Δx` | `(N_oct,1,n_diag)` | x-streaming (broadcast over ng) |
| `s_y = 2|μ_y|/Δy` | `(N_oct,1,n_diag)` | y-streaming |
| `Q` | `(N_oct or 1, ng, n_diag)` | weight-normalised per-ordinate source (anisotropy already folded in — see §5) |
| `in_x`, `in_y` | `(N_oct, ng, n_diag)` | incoming x-/y-face flux |
| `ψ̄ = psi_avg` | `(N_oct, ng, n_diag)` | cell-average flux (the SI/Krylov iterate datum) |
| `out_x = 2ψ̄ − in_x` | `(N_oct, ng, n_diag)` | outgoing x-face flux |
| `out_y = 2ψ̄ − in_y` | `(N_oct, ng, n_diag)` | outgoing y-face flux |

### 1.2 The reframe — x-scan, y-march

Fix a y-row `j`. March down the octant's x-direction, treating `in_y = ψy_in(i,j)`
as a GIVEN (it was produced by the previous y-row's `out_y`). Then the x-face
recurrence is, substituting `out_x = 2ψ̄ − in_x` with `ψ̄ = (Q + s_x·in_x +
s_y·ψy_in)/D` and `in_x(i+1) = out_x(i)`:

```
out_x = 2·(Q + s_x·in_x + s_y·ψy_in)/D − in_x
      = (2·s_x/D − 1)·in_x  +  2·(Q + s_y·ψy_in)/D
      = α·in_x + β
```

so the SCAN COEFFICIENTS are (VERIFIED numerically at the last-ULP =
principled-equiv):

> **α(i) = 2·s_x(i)/D(i) − 1**,   **D = σ_t + s_x + s_y**
> **β(i) = 2·(Q(i) + s_y(i)·ψy_in(i)) / D(i)**

and the per-row x-face sequence is `ordinate_scan(α[:,row], β[:,row], ψ_x_in_bc)`
(`orpheus/sn/spatial/scan.py:80-185`) — the EXACT same `ψ[i+1]=α[i]ψ[i]+β[i]`
Blelloch closed form the 1-D path calls. The cell average + transverse out-face
are then recovered PER CELL from the produced faces:

> **ψ̄(i) = 0.5·(in_x(i) + out_x(i))**   — see §5 anti-surprise (1): this STILL
> holds in 2-D, identically (it is `0.5·(in_x + (2ψ̄−in_x)) = ψ̄`); the transverse
> coupling is ENTIRELY inside β, NOT in the averaging.
> **out_y(i) = 2·ψ̄(i) − ψy_in(i)**     — the y-out face threaded to the next row.

### 1.3 Degeneracy check (the unification proof)

At `s_y = 0` (no transverse axis, i.e. d=1): `D = σ_t + s_x`, `α = 2s_x/(σ_t+s_x) −
1`, `β = 2Q/D` — which is EXACTLY the slab cache `a_attenuation = 2s/(s+σ_t)−1` and
`b = 2·QV·inv_denom` (VERIFIED: `2-D alpha(sy=0) == slab a` bit-exact). So
**`CumprodScan` IS `ScanMarch` at d=1** (Fork A2 in the verification plan) — the
1-D `_sweep_1d_unified` slab joint-batch IS `scan(x)` with no march. This is the
Cardinal-Rule-2 concept-count payoff (#222 acceptance bullet 6).

### 1.4 Operation-order / bit-identity caveat

The crosswalk reproduces ψ̄ and the out-faces to **FP-association** (nulp), NOT
bit-identity vs the wavefront — confirmed: `out_x direct = 1.5489795918367344` vs
`α·in_x+β = ...346` (2 ULP). This is BY DESIGN (#222: row-march vs anti-diagonal
are different valid linearizations of the same triangular solve) and is the
explicit "do NOT demand `array_equal` across schedules" boundary
(`scan_march_verification.md` G2, `_NULP_BOUND=128`). The kernel's own left-fold
`((σ_t+s_x)+s_y)` stays bit-id-load-bearing WITHIN a schedule.

---

## 2. `a_attenuation` CACHE STATUS + THE 2-D EXTENSION

### 2.1 Today: 1-D ONLY

`CollisionCache` (`orpheus/sn/spatial/sweep_cache.py:337-458`):
- `inverse_denom`, `a_attenuation`, `cumprod_a` — each `(N, ng, nx)`, chain-ordered
  along the cell axis (axis 2), built by `from_geometry(geom, sig_t)` (`:388-458`).
- **`a_attenuation = 2|μ|·A_total · inverse_denom − 1`** (`:446-449`); for slab
  `A_total=2, A_down=1, dA_w=0` ⟹ `a = 2s/(s+σ_t)−1` (the Cartesian form #222 Q2
  cites: `a = (s−σ_t)/(s+σ_t)`, identical).
- Built/memoised on the MESH: `_ensure_coll_cache` stashes it as `sn_mesh._coll_cache`
  (`sweep.py:387-410`), `GeometryCoefficients` as `sn_mesh._geom_cache`
  (`sweep.py:378-384`); both built once per `(mesh, σ_t)` epoch (the
  `_build_count==1` invariance pin, `sweep_cache.py:371-377`).
- **Consumed ONLY by `_run_1d_sweep`** (`sweep.py:587-588`, `:737-739`): `coll.a_attenuation[ords]`
  → `ordinate_scan`. The `GeometryCoefficients.from_mesh_and_quad` populator ASSERTS
  `reduced is not None` (`sweep_cache.py:217-221`) — it literally REFUSES to build for
  a 2-D Cartesian mesh ("2-D Cartesian wavefront uses anti-diagonal scheduling, not
  the chain scan"). So the cache is structurally 1-D today.

### 2.2 The 2-D wavefront does NOT use it (the #206 inline recompute)

`SweepDependencyGraph.apply_windowed` (`sweep_graph.py:849-927`) computes `D` and the
closure INLINE every level via `cell_update.cell_kernel_batch(...)` (`:905-916`) — `denom
= σ_t + s_x + s_y` rebuilt per cell, per level, per sweep, every source iteration.
The matvec twin `residual_windowed` (`:929-969`) does the same via
`residual_kernel_batch`. This is precisely the #206 duplication ("matvec/sweep both
recompute `denom = 2|μ|·A_down + dA_w·c_out + σ_t·V`"). #222 Q2 / `scan_march_verification.md`
decision (c) target: cache the flux-independent `D`/`α` so BOTH paths read one source.

### 2.3 The 2-D extension (decision (c) — per-line, mesh-memoised, `(d−1)`-slab)

The 1-D cache key is `(N, ng, nx)` because the chain is 1-D. At 2-D the denominator
`D[n, g, i, j] = σ_t[g,i,j] + s_x[n,i] + s_y[n,j]` depends on BOTH axes, so a naive
full cache is `(N, ng, nx, ny)` — which would BLOW the `(d−1)`-slab peak-memory win
(`scan_march_verification.md` G7.b explicitly catches a "materialize full-grid D"
regression). Two structural facts make a per-LINE cache the right shape:
- **`α` is flux-INDEPENDENT and iteration-INVARIANT** (depends only on `s_x, s_y, σ_t`
  — all fixed within a `(mesh, σ_t)` epoch). So it can be memoised once, like the 1-D
  `a_attenuation`.
- **The scan is PER-LINE** (`scan(x)` for a fixed `j`): the working set the scan needs
  is `α[:, :, :, j]` / one y-row — `(N, ng, nx)`, the `(d−1=1)`-slab. β IS flux-dependent
  (it carries `s_y·ψy_in`), so β stays per-iteration on the SWEEP (the #222 Q2 placement
  table: `b` → sweep), recomputed per line per row.

**Recommended 2-D cache shape:** `α[n, g, i, j]` and `inv_denom[n, g, i, j]` are the
flux-independent quantities; to keep peak memory `O((d−1)-slab)` per sweep, either
(i) mesh-memoise the full `(N, ng, nx, ny)` `α`/`inv_denom` ONCE per `(mesh, σ_t)`
(amortised across all iterations — the memory is `2× (N·ng·nx·ny)` static, paid once,
NOT per-sweep transient; the per-SWEEP working set stays the per-line `(N,ng,nx)` slice
+ the rolling β), OR (ii) recompute `α`/`inv_denom` per LINE inside the march (cheaper
memory, more flops). The verification plan's G7.b pins the per-SWEEP working set as
`(d−1)`-slab regardless; the static mesh cache is a `FOUNDATION` Pattern-2 pin (the
cached 2-D `α` term-equals the inline `cell_kernel_batch` `2s_x/D−1`). **Subsumes the
#206 inline recompute** — make the cached `α` the single source for sweep AND matvec.

⚠ Build-site nuance: the 1-D cache is populated by walking `dag_walk` (slow, once);
the 2-D `α` is a pure broadcast `2·s_x[:,:,None] / (σ_t + s_x[:,:,None] + s_y[:,None,:])
− 1` over `streaming_x`/`streaming_y`/`σ_t` (no chain walk needed — Cartesian has a
trivial cell order). So a NEW `from_cartesian`-style populator, NOT the `from_mesh_and_quad`
one (which asserts `reduced is not None`).

---

## 3. MARCH + BOUNDARY-SHED SEAM

### 3.1 How y-rows thread ψy (the march)

The y-coupling is `out_y(i,j) = 2·ψ̄(i,j) − ψy_in(i,j)`, and `ψy_in(i,j+1) =
out_y(i,j)` (the downstream y-face of row `j` is the upstream y-face of row `j+1`),
threaded in the octant's y-sweep direction (`sy>=0` → ascending `j`, `sy<0`
descending — the `face_in_y`/`face_out_y` offsets, `sweep_graph.py:504-505`). Row
`j=0`'s `ψy_in` is the domain y-inflow boundary trace (`boundary.face_view(y_in_face)`).
So the march is: for each y-row in octant-direction, (a) build β with the current
`ψy_in` row, (b) `ordinate_scan` the x-faces, (c) recover ψ̄ + `out_y`, (d) carry
`out_y` as the next row's `ψy_in`. This is `scan(x) ∘ march(y)`.

The **per-octant independence** is structural: the wavefront has one
`SweepDependencyGraph` per in-plane octant (`OctantLabel`, 4 in 2-D), each an
independent DAG; the scan-march loops the same 4 octants (the outer loop in
`sweep_octant_group`, `sweep.py:874-957`, becomes "for octant: scan-march"). No
cross-octant data dependence in the bare sweep (the boundary coupling is the sibling
`−B`, externalised — §3.3).

### 3.2 Reflective-BC outflow shedding — the seam vs the wavefront

The wavefront sheds domain-edge outflow PER LEVEL, ι*-style: `_MovingFrontier.shed`
(`sweep_graph.py:311-326`) captures the high-edge cells' out-faces into per-axis
`capture_x`/`capture_y` as the frontier advances; `sweep_octant_group` then writes
them into the trace (`sweep.py:956-957`):

```
boundary_flux.face_view(x_out_face)[oct_idx] = capture_x   # (N_oct, ng, ny)
boundary_flux.face_view(y_out_face)[oct_idx] = capture_y   # (N_oct, ng, nx)
```

For the **scan-march** the shed is simpler and direct (no rolling frontier):
- **x-outflow** = the LAST scan value of EACH row's x-chain = `out_x[-1, j]` for every
  `j` ⟹ the full `(N_oct, ng, ny)` x-outflow face is the stack of per-row last-scan
  values (one per y-row). [If `sx<0`, "last" is the low-x edge — octant-direction.]
- **y-outflow** = the LAST row's `out_y` over all `i` = `out_y[i, j_last]` ⟹ the full
  `(N_oct, ng, nx)` y-outflow face is the final row's transverse out-faces. [If `sy<0`,
  the first row processed.]

Both write into `boundary_flux.face_view(<out_face>)[oct_idx]` exactly as the
wavefront does — the SHED TARGET is identical, only the SOURCE (last-scan / last-row
vs per-level capture) differs.

### 3.3 ERR-056 shared-face-shed-order concern

ERR-056 (the diagonal-cubature shared-face reflect order, Phase 3 / `project_wave_o`)
is the risk that when a face is touched by >1 octant, the reflect must happen after
the LAST outflowing group. In the BARE sweep the shed write is disjoint-by-ordinate
(distinct octants own disjoint ordinate slices of a face — `sweep.py:950-957` notes
the write "touches only this octant's ordinates"), so the scan-march's shed is
ι*-faithful the same way. The ORDER concern lives in the Gauss-Seidel resolvent's
`reflect` step (`_sweep_2d_scheduled`, `sweep.py:1062-1068`), which is SCHEDULE-level,
NOT sweep-internal — the scan-march inherits it unchanged (it is a different x/y
SCHEDULE inside the bare sweep, orthogonal to the octant-group reflect schedule).
The verification plan's G4.b pins this; flags it is only fully exercised at d≥3 in
pure 2-D Cartesian.

⚠ Note: the bare sweep does NOT re-apply the reflective `R·G` (Wave O O.4b
BC-extraction — the sweep reads `ψ.boundary.inflow` as a GIVEN and writes
`ψ.boundary.outflow`; the reflection is the sibling `B` operator). So the scan-march
shed writes the OUTFLOW; the inflow seed is read from the trace `boundary.face_view`.
This matches `apply_windowed`'s contract exactly — no new BC logic.

---

## 4. THE STRATEGY SLOT

### 4.1 Where `ScanMarch` slots

`orpheus/sn/sweep_strategy.py`. The Protocol (`:164-217`) requires `sweep`,
`residual`, `residual_transpose`, `supports`. `ScanMarch` subclasses `_SweepStrategy`
(`:225`, the frozen-dataclass base with the `__post_init__` compatibility guard).
Registry: `SWEEP_STRATEGIES` tuple (`:488-492`) — `default_for` (`:495-519`) returns
the FIRST whose `supports(mesh).ok`. Selection priority decides Fork B1 (opt-in:
`ScanMarch` AFTER `MovingFrontierWindow`) vs Fork B2 (default: `ScanMarch` FIRST).

### 4.2 What `CumprodScan` does today (the d=1 unification target)

`CumprodScan` (`:287-334`): `supports` = `Compatibility(mesh.is_1d, ...)`; `sweep`
delegates to `_sweep_1d_unified(Q, sig_t, mesh, bf, initial_guess=...)` (`:320-322`)
with a `moment_projection` REFUSAL guard (1-D can't do moment output — the M-M Carlson
seed reads per-ordinate, L21, `:310-319`); `residual` → `operator._apply_1d` (`:328`);
`residual_transpose` → `operator._apply_1d_transpose` (`:334`). Under Fork A2 (#222
unification), `ScanMarch` ABSORBS this: `ScanMarch(d=1).sweep` IS the `scan(x)` form
(= `_sweep_1d_unified`'s slab joint-batch / the curvilinear per-ordinate scan), so
`CumprodScan` retires INTO `ScanMarch` (test migration: every `CumprodScan(mesh)` ref
rewires — `test_wavefront_cumprod_equivalence.py`, the registry, `default_for`; zero
dangling refs is the audit deliverable). ⚠ CURVILINEAR (sphere/cyl) stays a d=1
`ScanMarch` concern with the M-M angular thread folded into β (it ALREADY is — the 1-D
curvilinear `b = 2·(QV + ang_contrib)·inv_denom`, `sweep.py:739`); the 2-D
`scan(x)∘march(y)` is Cartesian-only (no curvilinear 2-D — the anti-hyperplane lattice
is a Cartesian object, `geometry.py:677-682`).

### 4.3 The `supports` predicate

`ScanMarch` is the UNION of the cases `CumprodScan` (1-D any geometry) and the 2-D
window (Cartesian d≥2) cover — but the two have DIFFERENT predicates:
- 1-D curvilinear (sphere/cyl): `is_1d` AND NOT `is_cartesian` — the chain scan with
  the M-M β. Valid.
- 1-D Cartesian (slab): `is_1d` AND `is_cartesian` — `scan(x)`, no march. Valid.
- 2-D Cartesian: `is_cartesian` AND `ndim==2` — `scan(x)∘march(y)`. Valid.
- 1-D NON-Cartesian 2-D? does not exist (curvilinear is 1-D only).

So the correct predicate is **`mesh.is_1d OR (mesh.is_cartesian and mesh.ndim >= 2)`**
— equivalently `mesh.is_1d or mesh.is_cartesian` (since `is_1d` already covers d=1 of
either geometry, and `is_cartesian` covers Cartesian d≥1; the union is "1-D anything,
OR Cartesian anything"). ⚠ The brief's proposed "`is_1d OR is_cartesian`" is CORRECT
and is the cleanest form: it admits 1-D-any-geometry (via `is_1d`) AND Cartesian-any-d
(via `is_cartesian`), and EXCLUDES the only impossible case (curvilinear d≥2, which is
not constructible anyway). `is_cartesian` (`geometry.py:643-658`) = `curvature is None`;
`is_1d` (`:628-641`) = `ndim==1`. The two are ORTHOGONAL (a slab is both; a cylinder is
`is_1d` not `is_cartesian`; a 2-D Cartesian is `is_cartesian` not `is_1d`) — confirmed
by the S0 memo's `reduced ⟺ is_1d` coincidence finding and the orthogonality pin.

⚠ d≥3 Cartesian: `is_cartesian` admits it, so `ScanMarch.supports` would return ok —
matching the recursive `scan(x)∘march(y,z)` design (#222 scope), verified synthetically
(G7.a) but with the d=3 contiguity/perf deferred. `MovingFrontierWindow` stays `d==2`
in `supports` (its orchestrator is 2-D until C3.6); `ScanMarch` can be the genuine
d-general production primitive alongside `FullFieldWavefront` (the oracle).

### 4.4 The matvec twin (`residual`)

L21 demands the apply-direction scan-march mirror the sweep. The 1-D twin is
`operator._apply_1d` (`operator.py:1464-1493`); the 2-D twin is `_apply_2d_cartesian`
(`:1641-1789`, routing `graph.residual_windowed`). The scan-march `residual` is the
apply-direction recurrence: `residual_kernel_batch` (`diamond.py:331-365`) computes
`r = D·ψ̄_probe − (Q + s_x·in_x + s_y·in_y)` and `out_a = 2·ψ̄_probe − in_a` from the
PROBE — i.e. the SAME α/β line structure but evaluating the residual at a given probe
ψ̄ rather than solving. ⚠ A2D-1 source-hash pin (`scan_march_verification.md` L3):
add the scan-march matvec as a NEW method (e.g. `_apply_2d_cartesian_scanmarch`),
leaving `_apply_2d_cartesian` UNTOUCHED → the hash stays free-green (Fork B1). The
`residual_transpose` (adjoint) for 2-D Cartesian is DEFERRED (raises
`NotImplementedError`, `sweep_strategy.py:364-378`) — the scan-march must preserve that
deferral boundary (an unsupported adjoint is a deferral raise, never a silent wrong
answer).

---

## 5. ANTI-SURPRISES (what would make the simple α·ψ+β recurrence WRONG)

Everything the scan-march must ALSO handle to be a faithful `FullFieldWavefront`
replacement at d=2:

**(1) ⭐ `psi_avg = 0.5·(in_x + out_x)` STILL HOLDS in 2-D — NOT a surprise, but the
COMMON wrong assumption.** One might think the 2-D cell average needs the y-faces too.
It does NOT: from the WDD closure `out_x = 2ψ̄ − in_x`, `0.5·(in_x + out_x) = ψ̄`
IDENTICALLY (verified). The transverse y coupling enters ONLY through β (it shifts ψ̄
via the source), so once the x-scan produces the x-faces, ψ̄ is recovered by the SAME
1-D averaging. The y-out face `out_y = 2ψ̄ − ψy_in` is then computed from ψ̄. ✅ So the
scan-march genuinely reduces to per-row 1-D scans + a thread — no hidden 2-D coupling
in the averaging.

**(2) ⭐ Anisotropic moment coupling (P1+) is OUTSIDE the recurrence.** `Q` is the
per-ordinate `AngularSourceSink` — the scattering operator's `R·Λ·M` (moment path,
`scattering.py:161-239`) and fission ALREADY projected the anisotropic source to
per-ordinate magnitude BEFORE the sweep (`transport_sweep` consumes ONE
`AngularSourceSink`, `sweep.py:104-245`; the `/W` projection is producer-side, lesson
L18). So the recurrence is per-ordinate-SCALAR (in `(N_oct, ng)` lanes) — the P1+
anisotropy never enters the cell kernel; it is folded into `Q` upstream. ✅ The
scan-march handles the SAME Q the wavefront does. (The Mode-9 gate G4 runs P1-aniso to
catch a schedule that secretly changes the operator — but the recurrence itself is
aniso-blind.)

**(3) ⭐ Per-ordinate vs moment-windowed OUTPUT (Phase 5a/5c) is an OUTPUT projection,
NOT a recurrence change.** `apply_windowed` branches angular-write vs moment-`+=`
(`sweep_graph.py:919-927`): in moment mode it projects `ψ̄` per anti-diagonal into the
harmonic tensor (`moment_buf += einsum("nlm,ngd,n->lmgd", Y, ψ̄, w)`) INSTEAD of
materialising the angular field (the ~3× peak-mem win). The scan-march must replicate
this: after recovering `ψ̄` per row, project into `moment_buf` if `moment_projection`
given (per-row `einsum` instead of per-anti-diagonal). The RECURRENCE is unchanged;
only the OUTPUT accumulation differs. ⚠ The moment guard: `CumprodScan.sweep` REFUSES
moment output (1-D, `:310-319`); the 2-D scan-march MUST support it (it is the windowed-SI
peak-mem path). The cross-octant/cross-row `+=` reorders the ordinate sum ⟹
principled-equiv not bit-id (already the case for the wavefront, ≤4 ULP).

**(4) Multi-group coupling: there is NONE inside the recurrence.** The cell kernel is
group-diagonal (`σ_t[g]`, `Q[g]` — the `ng` axis is a passive broadcast lane;
`cell_kernel_batch` carries `ng` as a vectorised dimension, no inter-group term). Group
coupling is in the SCATTERING source (upstream, in Q) and the outer iteration, NOT the
sweep. ✅ The scan processes all `ng` as parallel lanes (`ordinate_scan`'s trailing-lane
contract, `scan.py:106-119`) — same as the 1-D path's `(nx, K, ng)` joint batch.

**(5) Octant / sign handling.** 4 in-plane octants (2-D), each an INDEPENDENT DAG /
scan-march (no cross-octant coupling in the bare sweep — the reflective coupling is the
sibling `−B`, §3.3). The x-scan direction = `sign_x` (ascending-i for `sx≥0`,
descending for `sx<0` — `face_in_x`/`face_out_x` offsets, `sweep_graph.py:504-505`); the
y-march direction = `sign_y`. The pure-z degenerate octant (`sx==sy==0`, 2-D only) is
`ψ̄ = Q/σ_t`, NO faces, NO scan (`sweep.py:887-902`) — the scan-march must keep this
short-circuit. ✅ Mirrors `sweep_octant_group`'s octant loop.

**(6) Curvilinear pole (`a=0`, `D` resonance) is d=1-ONLY.** The cylindrical pole
`a_attenuation[i]=0` reset (ERR-054/#209) and the M-M angular thread (β's `ang_contrib`)
are CURVILINEAR — i.e. `ScanMarch(d=1)` territory (the Fork-A2 absorbed `CumprodScan`
curvilinear path). The 2-D `scan(x)∘march(y)` is Cartesian-ONLY (no curvilinear 2-D),
so the pole reset NEVER occurs in the 2-D march. BUT: `ordinate_scan` already handles
`a=0` (and the ERR-057 gradual denormal underflow) via the finite-test dispatch
`if not np.all(np.isfinite(closed_form)): return _pair_monoid_scan(...)`
(`scan.py:176-185` — the ERR-057 fix is LANDED; the test-architect's "denormal gap" is
already closed in this worktree). So the scan-march inherits the robustness for FREE
by reusing `ordinate_scan` per line. ✅ No new underflow handling needed; the
conditioning gate G3 is inherited.

**(7) The `b`/source `2×` factor + cell-volume.** The 1-D path scales `QV = Q·V`
(`sweep.py:484`) then `b = 2·QV·inv_denom`. The 2-D kernel's `Q_cells` is ALREADY
weight-normalised AND the 2-D denominator has NO explicit `V` (it is `σ_t + s_x + s_y`,
the `Δ`-normalised form where `s_a = 2|μ|/Δ_a`). So the 2-D β = `2·(Q + s_y·ψy_in)/D`
uses the Δ-normalised `Q` directly (NO `·V`) — the `V` scaling is the CURVILINEAR/1-D
chain convention (`denom = 2|μ|·A_down + σ_t·V`), the Cartesian wavefront uses the
already-divided-by-Δ form. ⚠ Build-discipline: do NOT copy the 1-D `QV = Q·V` into the
2-D march — the 2-D `Q_cells` convention is the `cell_kernel_batch` one (`Q + Σ s_a·in_a`
over `σ_t + Σ s_a`), NOT the `2·QV·inv_denom` chain form. Get the convention from
`cell_kernel_batch` (`diamond.py:322-327`), not from `_run_1d_sweep`. The factor-of-2
DOES appear (β has the `2×`) but the V does not.

---

## KEY FILE:LINE INDEX (re-confirm via Nexus at pickup; SHAPE durable)

| seam | file:line |
|------|-----------|
| 2-D DD kernel (α/β ground truth) | `orpheus/sn/spatial/diamond.py:273-329` (`cell_kernel_batch`), `:331-365` (`residual_kernel_batch`) |
| `ordinate_scan` (the reused scan; ERR-057 finite-test guard LANDED) | `orpheus/sn/spatial/scan.py:80-185` (`:176` the `isfinite` dispatch) |
| `a_attenuation` cache (1-D only) | `orpheus/sn/spatial/sweep_cache.py:438-449`; `from_geometry` `:388-458`; `assert reduced is not None` `:217-221` |
| 1-D sweep (β-build = §1) | `orpheus/sn/sweep.py:413-779` (`:593` slab b, `:739` curvilinear b) |
| 2-D orchestrator (octant loop, capture shed) | `orpheus/sn/sweep.py:796-957` (`:956-957` the shed) |
| windowed walk (inline D recompute = #206) | `orpheus/sn/sweep_graph.py:849-927` (`apply_windowed`), `:929-969` (`residual_windowed`); `_MovingFrontier.shed` `:311-326` |
| `s_a = 2|μ|/Δ` accessor | `orpheus/sn/geometry.py:660-687` (`streaming(axis)`), `:1347-1354` (`streaming_x`/`y` build) |
| `is_1d` / `is_cartesian` | `orpheus/sn/geometry.py:628-641` / `:643-658` |
| strategy slot + registry + `default_for` | `orpheus/sn/sweep_strategy.py:287-334` (`CumprodScan`), `:488-519` (registry/`default_for`) |
| 2-D matvec twin (A2D-1 hash, WRAP not relocate) | `orpheus/sn/operator.py:1641-1789` (`_apply_2d_cartesian`); 1-D twin `:1464-1493` |
| V&V gate plan (companion) | `.claude/agent-memory/test-architect/scan_march_verification.md` |
