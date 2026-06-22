---
name: issue-196-phase-g-step2-5-further-collapse
description: Read-only feasibility study of four further-collapse options after Phase G Step 2.5 landed (11 -> 8). Verdicts per Q1-Q4 with Pattern 6 vs Pattern 2 spectrum scoring and a realistic concept-count floor.
metadata:
  type: project
---

# Phase G Step 2.5 further-collapse feasibility (Issue #196)

Read-only investigation of four candidate unifications beyond the
Step 2.5 landing (commits `23766ee..1a699b0`). Outcome: **two of
the four are achievable, two are principled deferrals**. Realistic
concept-count floor is **11 -> 6 (~1.8:1)**, not 11 -> 3.

Citations are file:line into `refactor/sn-operator-algebra`.

---

## Q1. Can `_sweep_1d_cartesian` and `_sweep_1d_curvilinear` share a skeleton?

### Side-by-side anatomy

| Stage | `_sweep_1d_cartesian` (sweep.py:171-289) | `_sweep_1d_curvilinear` (sweep.py:296-523) |
|---|---|---|
| Outer iteration | `for n in range(N)` (line 241) | `for p_idx, level in enumerate(levels)` + inner `for m_local, global_n` (lines 418, 443) |
| Carlson half-angle seed pass | **absent** | per-level call to `carlson_inward_sweep_from_source` (line 431) producing `psi_angle = phi_aux.T.copy()` |
| Spatial-upstream init | BC inflow at left/right face (lines 252-255) | Outward: pole-face IC `psi_bc[pole_key][global_n]` (line 461). Inward: BC inflow (line 454) |
| Inner fold | `for visit in iter_cell_visits(ordinate_idx=n)` (line 263) | `for visit in iter_cell_visits(ordinate_idx=ordinate_idx, mu_level_idx=...)` (line 483) |
| UpstreamState construction | `angular_upstream=None` (line 267) | `angular_upstream=psi_angle[i]` (line 490) |
| Per-cell update call | identical: `cell_update.update(visit, ..., upstream)` (lines 269-274 vs 492-497) | same |
| Angular-state thread | absent (slab `outgoing_angular_state is None`) | `psi_angle[i] = result.outgoing_angular_state` (line 501) |
| Spatial closure thread | unconditional `psi_spatial_in = result.outgoing_spatial_flux` (line 280) | guarded by `if result.outgoing_spatial_flux is not None` (line 507) |
| BC outflow capture | unconditional (lines 284-287) | guarded by `abs_mu_n >= eps` non-degenerate (lines 514-516) |
| Persistent-state cache | reflective face buffers only | + pole_key + phi_prev_key (lines 520-521) |

### Shared structure is real

The inner fold body is **already structurally identical** at the call
site granularity. The DD strategy at `diamond.py:128-175` consumes both
paths in one body — that part of the collapse already shipped. The
visible differences are at **outer-loop level** (per-level vs
per-ordinate iteration) and at **pre/post passes** (Carlson seed
generation + state caching).

### Strategy-injection skeleton — feasibility

`PsiHalfAngleSeed` and `PoleAngularClosure` Protocols exist
(`psi_half_angle_seed.py`, `pole_angular_closure.py`). **However, they
are matvec-side Protocols, not sweep-side stages.** Evidence:

- `PoleAngularClosure.__call__` consumes `(psi_cells, alpha_half,
  redist_dAw, tau_mm, volume, level_indices, carlson_context)` and
  returns a **redistribution tensor** (`pole_angular_closure.py:266-276`).
  That is the operator-form L·psi computation. It is **not** what the
  sweep does — the sweep's M-M redistribution is computed implicitly via
  `psi_angle[i]` thread inside `DiamondDifference.update` and the
  WDD spatial closure.
- `PsiHalfAngleSeed.__call__` takes `(psi_level, context)` and returns
  `(ng, nx)` — that IS what the sweep's pre-pass needs, but the sweep
  doesn't currently invoke it. It calls
  `carlson_inward_sweep_from_source` directly (sweep.py:431). The
  reason is that the sweep generates `Q_bar` from
  `phi_0_prev` (the previous-iteration scalar flux cached in
  `psi_bc[phi_prev_key]`), not from an angular-flux level tensor
  `psi_level` — there is no `psi_level` at sweep time, only a Carlson
  pre-pass that bootstraps the angular state for the level.

The asymmetry is structural: matvec has `psi` in hand and projects
moments; SI has `phi_0_prev` cached and skips the projection. They
both call `carlson_inward_sweep_from_source` — that's the single
source of truth (Pattern 2 — already achieved).

### Verdict

**Partial collapse is achievable, but the "strategy-injection
skeleton" framing is the wrong abstraction.** The genuinely-different
stages are:

1. **Level enumeration** — `[None]` (slab/sphere) vs
   `quad.level_indices` (cylinder). Data-driven (Pattern 2). Already
   uniform within `_sweep_1d_curvilinear`. Slab can adopt the same
   shape with `levels = [None]`.
2. **Carlson seed pass** — present for curvilinear, absent for slab.
   The natural unification: factor the seed pass into a free function
   `_seed_angular_state(reduced, phi_0_prev, ...) -> ndarray | None`
   that returns `None` for slab. The sweep then has a uniform
   `psi_angle = _seed_angular_state(...)` call where slab gets `None`
   and constructs `UpstreamState(angular_upstream=None)` from it.
3. **Pole-face IC** vs **BC-only IC** — both already encoded as
   `psi_spatial_in = <something>(direction_sign, ordinate)`. Trivially
   data-driven.
4. **Persistent caches** — `psi_bc[pole_key]`, `psi_bc[phi_prev_key]`
   are needed only when `requires_upstream_angular_state`. The cache
   dict is already keyed by string; one guarded `if reduced.requires_upstream_angular_state`
   suffices.

**Sketch of the unified body (~30 lines):**

```python
def _sweep_1d_unified(Q, sig_t, sn_mesh, psi_bc, Q_aniso):
    reduced = sn_mesh.reduced
    has_angular_dag = reduced.requires_upstream_angular_state
    quad = sn_mesh.quad
    # Level enumeration — sphere=[None], cylinder=level_indices,
    # slab=[None] (no per-level angular DAG; psi_angle stays None)
    levels, level_ord_lists = _enumerate_levels(quad, reduced)
    inflow_full = _resolve_outer_inflow(sn_mesh, psi_bc)
    angular_flux, scalar_flux = _alloc_outputs(sn_mesh, Q, quad)
    for p_idx, level in enumerate(levels):
        psi_angle = (_seed_carlson(reduced, sn_mesh, psi_bc, level)
                     if has_angular_dag else None)
        for global_n in level_ord_lists[p_idx]:
            psi_spatial_in = _resolve_inflow(sn_mesh, psi_bc, global_n,
                                              inflow_full)
            for visit in sn_mesh.iter_cell_visits(...):
                up = UpstreamState(
                    spatial_upstream=psi_spatial_in,
                    angular_upstream=(psi_angle[visit.cell_idx]
                                       if psi_angle is not None
                                       else None),
                )
                result = sn_mesh.cell_update.update(visit, ...)
                if psi_angle is not None and result.outgoing_angular_state is not None:
                    psi_angle[visit.cell_idx] = result.outgoing_angular_state
                if result.outgoing_spatial_flux is not None:
                    psi_spatial_in = result.outgoing_spatial_flux
                # accumulate angular_flux, scalar_flux
            _capture_outflow(sn_mesh, psi_bc, global_n, psi_spatial_in)
    if has_angular_dag:
        _cache_state(psi_bc, angular_flux, scalar_flux)
    return angular_flux, scalar_flux
```

Concept count: **2 -> 1** for the 1-D skeleton (the cartesian/curvilinear
fork at sweep.py:158-163 collapses to data-driven enumeration). 4-5
helper functions are introduced but they are named primitives — each
testable in isolation. The body reads `for level: for ordinate: for
visit:` — that's the math. Net win: ~290 LOC -> ~110 LOC, three
distinct sweep entry-points -> two.

**Pattern axis**: solid Pattern 2 win (single source of truth for the
fold). The risk is that slab now carries the level-enumeration scaffold
even though slab has only one trivial level — but the scaffold is
no-cost data (`levels = [None]`) and reading `for level in [None]: for
n in range(N): ...` is barely heavier than `for n in range(N): ...`.

---

## Q2. Can `_sweep_2d_wavefront` join the unified 1D skeleton?

### Structural difference

`_sweep_2d_wavefront` (sweep.py:540-708) iterates **topological
levels** (anti-diagonals via `SweepDependencyGraph.levels`), each
level passed as a `SweepCellSlice` to `update_batch`. The work unit is
a slice of `n_diag` cells x `N_oct` ordinates simultaneously.

1-D iterates **single cells**, threading scalar state (`psi_face_in`)
across cells.

### Fold-shape arithmetic

The two are **algebraically the same fold** — forward substitution on
block-triangular L. But the **work-unit cardinality** differs by
factor `N_oct * n_diag` (typically 4 x sqrt(nx*ny) ~ 20-80x at
production grids).

Concretely:

- 1-D inner fold body (sweep.py:263-280): one `update` call per cell;
  state is one `(ng,)` scalar threaded along
- 2-D inner fold body (sweep_graph.py:326-352): one `update_batch`
  call per level; state is `(N_oct, nx+1, ny, ng)` face buffers
  mutated in place by the batched closure

### Forcing 1-D into length-1 slices: cost estimate

For nx=80 (typical), ng=2, N=8 (S8 GL):

- Current 1-D: 80 cells x 8 ordinates = 640 per-cell `update` calls
  per sweep. Each builds an `UpstreamState`, a `CellResult` (3-field
  dataclass), and an `np.ndarray` of shape (ng,).
- 1-D-via-slice variant: 80 length-1 `SweepCellSlice` allocations per
  ordinate per sweep, each carrying 12 fields including index arrays of
  size 1. Inside `update_batch`, every numpy op is on shape
  `(1, 1, ng)` — fully dominated by Python/numpy dispatch overhead.

`SweepCellSlice` has 12 fields (cell_update.py:268-279); naive
allocation is ~10x heavier than `CellVisit` + `UpstreamState`. The
numpy-vectorisation benefit is **negative** at `n_diag=1` because the
Python dispatch wins by O(constant) per call. Crude estimate: ~3-5x
slowdown on 1-D production sweeps.

### What `dag_walk` could yield

Three options were considered:

- **(a) WorkUnit Union** — Iterator yields `CellVisit | SweepCellSlice`,
  consumer type-switches. Adds dispatch, reads worse than current
  twin code. Anti-pattern (stringly-typed dispatch).
- **(b) Always SweepCellSlice** — described above; ~3-5x slowdown for
  1-D, no gain.
- **(c) Polymorphic per-mesh-type** — `dag_walk` returns a different
  primitive per mesh, with two consumers (per-cell folder, per-level
  folder). Doesn't unify; documents the split.

### Verdict

**Pattern 6 deferral.** The 1-D fold over per-cell `update` and the
2-D fold over per-level `update_batch` are **algebraically the same**
but the work-unit shapes differ by O(n_diag) — forcing them into one
shape either pessimises 1-D (length-1 slices) or pessimises 2-D
(per-cell calls). The right abstraction here is the **DD strategy
itself**: `DiamondDifference` already implements BOTH `update` (per-cell)
and `update_batch` (per-level), and they are bit-identical at slot
granularity. That IS the unification — the algebra is shared via the
strategy; the iteration shapes legitimately differ.

The naming `update` vs `update_batch` is the current concept-count
cost. See Q4 for whether that pair can be collapsed.

**Pattern axis**: this is the canonical Pattern 6 "two working
implementations earn the right to be unified, but ONLY when the
unification is honest." Here the unification path actively
pessimises one side. Premature abstraction.

---

## Q3. Can `dag_walk` be the canonical primitive (not an alias)?

### Call-site map

`iter_cells_by_direction` consumers (all matvec):

- `operator.py:752, 753` (spherical, +1 and -1)
- `operator.py:1009, 1063` (cylindrical, +1 and -1 per level)

`iter_cell_visits` consumers (all SI sweep):

- `sweep.py:263` (slab 1-D Cartesian)
- `sweep.py:483-486` (curvilinear)

`dag_walk` (currently a thin alias at `geometry.py:425-466`): **zero
production consumers**; only doc-string examples.

### Why SI needs per-ordinate, matvec needs per-direction

Read `_iter_cartesian_visits` at `geometry.py:573-598` and
`_iter_spherical_visits` at `geometry.py:600-629`: the StreamingTerms
on the emitted `CellVisit` are built per `(cell_idx, ordinate_idx)`
via `self.reduced.streaming_terms(cell_idx=i, direction_idx=ordinate_idx)`.
For sphere, `direction_idx` selects the global ordinate; for
cylindrical, the within-level azimuthal index. The DD strategy reads
`visit.streaming_terms.abs_mu`, `mu`, `tau_mm`, `alpha_in/out`, etc.
— all per-ordinate quantities.

The matvec at `operator.py:752` consumes ONLY `visit.cell_idx` (line
788 `i = visit.cell_idx`; everything else is recomputed from `A[i+1]`,
`A[i]`, `V[i]`, `mu[outgoing_mask]`). It does NOT read
`visit.streaming_terms`. So feeding it a representative-ordinate
visit is safe.

### Could `dag_walk(ordinate_idx)` subsume both?

Yes — that's literally what `iter_cell_visits(ordinate_idx)` already
does. The current `dag_walk(direction_sign, mu_level_idx)` picks a
representative ordinate internally (geometry.py:746-756) and delegates
to `iter_cell_visits`. The alias relationship is:

- `dag_walk(direction_sign, mu_level_idx)` =
  `iter_cell_visits(representative_ordinate_for(direction_sign),
                     mu_level_idx)`

### Migration cost

A canonical `dag_walk(*, ordinate_idx=None, direction_sign=None,
mu_level_idx=None)` with the constraint "exactly one of
ordinate_idx, direction_sign is required" subsumes both. Then
`iter_cell_visits` and `iter_cells_by_direction` become deprecated
aliases.

Code touched: 4 call sites in operator.py + 2 in sweep.py = 6 sites.
Each is a 1-line `for visit in sn_mesh.dag_walk(...)` rewrite. Plus
the alias bodies in geometry.py (lines 463-465, 631-756) can be
deleted, saving ~125 LOC in geometry.py.

### Verdict

**Canonicalise.** `dag_walk` is currently a stub doing PR-friendly
naming work; making it canonical removes 2 redundant iteration
methods. This is a clean Pattern 2 win: one named primitive, two
parameterisations.

**Pattern axis**: Pattern 2 single-source-of-truth, very low risk —
the alias already documents the relationship. Migration is
mechanical.

Concept count: **3 -> 1** for iteration primitives. The closeout
audit reads "three iteration methods" — this one stroke removes
the appearance of three concepts; what was hidden behind the
"thin alias" wording becomes architecturally legible.

---

## Q4. Can `update` and `update_batch` collapse to one method?

### Algebra orientations

`update` (diamond.py:128-175): solves `denom * psi_avg = source +
numer_upstream` where `numer_upstream` accumulates `(2*abs_mu*A_down)
* psi_spatial_in + (dAw * alpha_in) * psi_angle_in`. Threading is
**scalar via UpstreamState**.

`update_batch` (diamond.py:214-292): solves `denom = sigt + sx + sy;
psi_avg = (Q + sx*psi_in_x + sy*psi_in_y) / denom` and **scatters**
`2*psi_avg - psi_in_x` into `psi_x` at outgoing-x-face indices via
`s.psi_x[:, s.face_out_x_idx, jj, :] = ...`.

Two structural differences:

1. **State threading mechanism.** `update` threads scalar
   upstream-state THROUGH RETURN VALUE (`outgoing_spatial_flux`
   captured by caller and rebound to `psi_spatial_in`).
   `update_batch` threads array face fluxes through IN-PLACE MUTATION
   of `psi_x/psi_y` at indices supplied by `SweepCellSlice`.
2. **Algebra form.** 1-D uses the `cell_balance_terms` helper which
   returns the geometry-blind `(denom, numer_upstream)` pair (works
   for slab/sphere/cylinder via neutral curvature). 2-D uses
   `denom = sigt + sx + sy; numerator = Q + sx*psi_in_x + sy*psi_in_y`
   — bit-identity-pinned to legacy sweep.py:858-867 (diamond.py:243-251
   docstring).

### Re-expressing `update_batch` as a fold over per-cell `update`

Mechanically possible: for each `(ii, jj)` pair on the level,
build a synthetic 2-D CellVisit + UpstreamState, call `update`,
scatter the output. But:

- The 2-D algebra form is `(Q + sx*psi_in_x + sy*psi_in_y) / (sigt
  + sx + sy)`. The 1-D `cell_balance_terms` form expects `numer_upstream
  = 2*abs_mu*A_down*psi_spatial_in + dAw*alpha_in*psi_angle_in` and a
  single `psi_spatial_in`. Forcing 2-D into this shape needs **two**
  spatial-upstream channels (x and y), which the current `UpstreamState`
  doesn't carry.
- Even with that extension, per-cell looping over n_diag cells at
  numpy-shape `(N_oct, ng)` LOSES the simultaneous vectorisation across
  `n_diag` cells that makes `update_batch` ~10x faster than the legacy
  per-ordinate loop. Concrete: nx=ny=80, S8 -> 79 levels avg 40 cells
  each, 4 octants. Current `update_batch`: 4*79 = 316 calls per sweep
  with `(8, 40, ng)` shape ops. Per-cell-fold: 4*79*40 = 12640 calls
  with `(8, ng)` ops — 40x more Python dispatch.
- The bit-identity contract on legacy sweep.py:847-871 (diamond.py:243-252)
  is enforced by tests; rewriting the operation order breaks it.

### Building the right primitive (Pattern 5)

Could a single method take a "cell or slice" argument? Concretely a
`CellOrSlice` union — DD dispatches internally on type. That just
moves the if-statement inside the strategy (anti-pattern: stringly-typed
dispatch).

A more principled framing: **the 1-D fold is `update_batch` on a
"slice of length 1 with a 1-D `SweepCellSlice` variant"**. But:

- 1-D's state is genuinely scalar-threaded — every cell depends on
  the previous one's output. 2-D's anti-diagonal levels are
  **mutually independent**. The fold structures differ in what state
  flows between iterations: 1-D has accumulator-fold (state threaded);
  2-D has parallel-batch-fold (state in persistent buffers, no
  level-to-level threading in the inner step).

That's not just shape — that's a real difference in the fold's
operator algebra. Fold-with-accumulator vs fold-with-parallel-reduce
are different patterns.

### Verdict

**Keep dual methods.** `update` and `update_batch` describe the same
math discretised across different work-unit cardinalities. The
single-cell variant must exist for 1-D's scalar-threaded fold; the
slice variant must exist for 2-D's anti-diagonal vectorisation; and
the bit-identity contract on the slice variant pins it to the
legacy operation order. The right Pattern 2 win was already
achieved: the **per-cell algebra** (cell_balance_terms helper) is
the single source of truth. The `update_batch` body is a separate
numpy form because the 2-D balance equation has a different surface
shape (two streaming channels + scatter mutation) — it cannot be
expressed as length-1 instances of `update` without losing
vectorisation.

**Pattern axis**: Pattern 6 deferral with strong justification. The
two methods are NOT "twin paths to the same product" (anti-pattern)
because each is the **uniquely-best implementation** for its work-unit
shape. The shared structural truth — that they implement the same
math — is captured in the diamond.py module docstring and the
`is_linear` ClassVar consistency.

Possible follow-up: if/when LD or EC closures want batched 2-D, they
will need to follow the same dual-method pattern. That's also fine —
it's a strategy contract, not duplication. Worth flagging in
`docs/theory/discrete_ordinates.rst` as the architectural decision.

---

## Overall concept-count audit

Post-Step-2.5 baseline (per L12 audit): 11 -> 8.

If the Q1 and Q3 unifications ship:

| Concept | Before | After Q1+Q3 |
|---|---|---|
| Sweep bodies | 3 (cartesian, curvilinear, 2d-wavefront) | 2 (unified-1d, 2d-wavefront) |
| Iteration primitives | 3 (dag_walk alias, iter_cell_visits, iter_cells_by_direction) | 1 (canonical dag_walk) |
| DD geometry branches | 1 (already unified Step 2.5) | 1 |
| Cell-update methods | 2 (update, update_batch) | 2 (keep — Q4 verdict) |
| **Total** | **9 (8 + 1)** | **6** |

Realistic floor: **11 -> 6 (~1.8:1)**. Not the 11 -> 3 (~3.7:1)
original ambition, but a real win over the current 11 -> 8.

The 1.4:1 -> 1.8:1 improvement is honest. The 11 -> 3 target was
**incorrectly priced** in the original Phase G plan — Q2 (2-D fold
shape) and Q4 (per-cell vs slice methods) were under-analysed and
their non-unification is the principled conclusion, not residual
implementation debt.

---

## Pattern 6 vs Pattern 2 spectrum

| Question | Verdict | Pattern 2 (SSoT) signal | Pattern 6 (premature) signal |
|---|---|---|---|
| Q1 1-D unify | **Unify** | Inner fold body already identical; differences are pre/post passes that factor cleanly | Slab carries `levels=[None]` scaffold it doesn't need; minor cognitive cost |
| Q2 2-D join 1-D | **Defer** | Same algebraic fold; appealing symmetry | Length-1 slices pessimise 1-D ~3-5x; fold-with-accumulator vs fold-with-parallel-reduce ARE different patterns |
| Q3 dag_walk canonical | **Canonicalise** | One named primitive, two parameterisations; eliminates 125 LOC of alias indirection | None — alias is already documented as a stub |
| Q4 update/update_batch | **Defer** | Both implement same math | Different state-threading mechanisms (return-value scalar vs in-place mutation); different fold patterns; bit-identity contract on the slice form |

Two of the four (Q1, Q3) are Pattern 2 wins — the unifications
**reduce concept count and increase clarity**. Two of the four (Q2,
Q4) are Pattern 6 deferrals — the unifications **multiply concept
count or pessimise performance** because the underlying algorithms
have genuinely different operator-algebra shapes.

The user's elegance signal (`[[feedback-elegance-causes-collapse]]`)
is satisfied at 11 -> 6: deeper unification crosses the line into
forcing genuinely different algorithms into one shape.

---

## Recommendation to main agent

Two further-collapse units fit a future Phase G Step 2.6:

1. **Unify 1-D sweep bodies behind level-enumeration data + factored
   Carlson seed helper.** Mechanical. ~290 LOC -> ~110 LOC.
   `transport_sweep` dispatch collapses to `if reduced is None:
   _sweep_2d_wavefront else _sweep_1d_unified`.
2. **Canonicalise `dag_walk`.** Take `*, ordinate_idx | direction_sign`.
   Migrate 6 call sites. Delete `iter_cell_visits` and
   `iter_cells_by_direction`. ~125 LOC saved in geometry.py.

Defer Q2 and Q4 with a Sphinx note in `docs/theory/discrete_ordinates.rst`
"Architectural deferrals" subsection explaining WHY (the algebraic
fold shapes genuinely differ; 11 -> 3 target was over-priced).
