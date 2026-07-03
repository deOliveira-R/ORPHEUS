# A3 — SN loss-operator "sweep machinery" fragmentation map

Scoping pre-read for the planned orientation×kernel unification carve. **Read-only**; no
code touched. `file:line` are current as of this read (2026-06-29, branch `main`) — treat
them as *re-derive-via-grep* markers; the **structural claims are the durable part**. Line
numbers in the original brief matched the current file, so the references below are live.

> **Nexus hazard honored:** the MCP active graph is a sibling worktree, so this map is built
> from Read+grep (the "why" lives in docstrings/comments Nexus doesn't index). No Nexus
> structural query was trusted.

---

## 0. The four faces and how they are REACHED (entry points)

The operator layer (`orpheus/sn/operators/streaming.py`) is the consumer; the
`LossRepresentation` protocol (`loss_representation/__init__.py:229`) is the seam.

| Face | Math | Reached from | Status |
|------|------|--------------|--------|
| `sweep` | `(L+C)⁻¹q` | `InvertibleOperator.solve` (streaming.py:761 → :1001) **and** the operator-free `transport_sweep` (__init__.py:1938 → `default_for(mesh).sweep`) | full coverage |
| `loss_action` | `(L+C)ψ` | `InvertibleOperator.apply` (streaming.py:712 → :738) | full coverage |
| `loss_action_transpose` | `(L+C)ᵀφ` | `InvertibleOperator.apply_transpose` (streaming.py:740 → :757) | **1-D only; multi-D Cartesian DEFERRED** |
| `sweep_transpose` | `(L+C)⁻ᵀb` | — | **DOES NOT EXIST** (grep across `orpheus/ tests/ docs/` = 0 hits; no `InvertibleOperator.solve_transpose` either) |

`streaming_action` / `streaming_action_transpose` (the σ-free `Lψ` / `Lᵀφ`, reached from the
bare `StreamingOperator.apply`/`apply_transpose`, streaming.py:423/:455) are **not** separate
faces — they are `loss_action(0, ·)` / `loss_action_transpose(0, ·)` single-sourced through the
base (`_LossRepresentation.streaming_action`, __init__.py:376/:383). The within-group matvec
is affine in σ, so σ=0 gives bare streaming.

**Consequence for the carve:** the unification targets a 2×2 grid {forward, reverse} ×
{solve, apply}, but **one cell is empty** (`sweep_transpose`) and **one is half-filled**
(`loss_action_transpose`: 1-D done, multi-D deferred). The unified walk must express a
*partial/deferred* grid, not a full symmetric 2×2.

---

## 1. THE COVERAGE MATRIX (the heart)

Registered representations (`LOSS_REPRESENTATIONS`, __init__.py:1879, in selection-priority
order): `CumprodScan`, `ScanMarch`, `MovingFrontierWindow`, `FullFieldWavefront`. The two
**shared walk frames** `_OctantWalk` (:717) and `_OneDimScanWalk` (:2143) are NOT registered —
they are the reusable traversal engines the representations delegate into.

| rep (file:line) | `sweep` (fwd solve) | `loss_action` (fwd apply) | `loss_action_transpose` (adj apply) | dims / curvature |
|---|---|---|---|---|
| **`_LossRepresentation`** base (:333) | `raise NotImplementedError` (:500) | `raise NotImplementedError` (:514) | `raise NotImplementedError` (:527) | — (never instantiated; supplies the construction guard `__post_init__` :539 + the σ-free `streaming_action` defaults) |
| **`CumprodScan`** (:976) | DELEGATES → `_OneDimScanWalk.sweep` (:1012); guards `moment_frame` (1-D ⇒ raise, moment output is 2-D only) | DELEGATES → `_OneDimScanWalk.loss_action` (:1028) | DELEGATES → `_OneDimScanWalk.loss_action_transpose` (:1042) — **IMPLEMENTED** | 1-D any geometry (slab/sphere/cyl), affine-scannable scheme (supports :987) |
| **`_DAGWavefront`** base (:1050) | abstract (subclass) | abstract (subclass) | **DEFERS** (:1091) `raise NotImplementedError("…multi-D Cartesian adjoint is deferred (O.2b lands the 1-D reverse sweep first)")` — inherited by BOTH subclasses | Cartesian d=2 (supports :1074) |
| **`MovingFrontierWindow`** (:1108) | `_sweep_jacobi(interior=self._sweep_interior)` (:1139) → `_sweep_scheduled` → **`_OctantWalk.sweep_group`**; interior (:1145) = `graph.walk_windowed × _CellSolve` (wraps `DiamondDifference.cell_kernel_batch`) | **`_OctantWalk.loss_action`** (:1226); interior (:1230) = `graph.walk_windowed × _CellResidual` (wraps `residual_kernel_batch`) | inherited **DEFERRED** | Cartesian d=2 (rolling (d-1)-frontier; walk is d=3-capable at graph layer but `supports` stays d=2). **Selectable peer**, not default since Fork-B2 |
| **`FullFieldWavefront`** (:1281) | `_sweep_jacobi(interior=self._sweep_interior)` (:1424); guards `moment_frame` (oracle has no moment mode); interior (:1430) = `graph.walk_full × _CellSolve` | **`_OctantWalk.loss_action`** (:1497); interior (:1501) = `graph.walk_full × _CellResidual` | inherited **DEFERRED** | Cartesian **ANY-d** (supports :1304 overrides to `is_cartesian`). The d-generic verification ORACLE / spine |
| **`ScanMarch`** (:1553) | 1-D → `_OneDimScanWalk.sweep` (:1650); 2-D → `_sweep_jacobi(interior=self._sweep_interior)` (:1658) → `_OctantWalk.sweep_group`; interior (:1664) = **hand row-march** `scan(x)∘march(y)` via `_scanmarch_row` + `scheme.cartesian_scan_coefficients` (NOT the DAG) | 1-D → `_OneDimScanWalk.loss_action` (:1781); 2-D → `_OctantWalk.loss_action` (:1782); interior (:1786) = row-march via `_x_scan_faces` + `reflect_scan_coefficients` + `residual_kernel_batch` | 1-D → `_OneDimScanWalk.loss_action_transpose` (:1857); 2-D → **DEFERRED** (:1858) `raise NotImplementedError` | 1-D any geometry + 2-D Cartesian **facewise-coupling** schemes only (DD/Step, NOT LD; supports :1599). **2-D Cartesian PRODUCTION DEFAULT** since Fork-B2 |
| **`_OctantWalk`** shared frame (:717) | `sweep_group` (:798) — SOLVE frame | `loss_action` (:858) — APPLY frame | **NONE** — no transpose method exists on this frame | multi-D Cartesian; both directions ride `_interior_walk` (:753) |
| **`_OneDimScanWalk`** shared frame (:2143) | `sweep` (:2160) → `_run` (:2816) — SOLVE = Blelloch prefix **scan** over chain order | `loss_action` (:2222) → `_apply_walk` (:2250) — APPLY = per-cell **loop** over `dag_walk_cell_indices`, forward | `loss_action_transpose` (:2602) — ADJOINT = per-cell loop over `dag_walk_cell_indices`, **`reversed`** + boundary swap + Carlson mirror + `angular_adjoint` | 1-D any geometry |

**One-line reading of the matrix:** `sweep` and `loss_action` are fully covered (1-D via
`_OneDimScanWalk`, multi-D via `_OctantWalk`). `loss_action_transpose` is **covered only in
1-D**; every multi-D path (`_DAGWavefront`, `ScanMarch` 2-D) raises a deferral. The base raises
"must implement". The 1-D `sweep` solve uses a **scan**, while the 1-D apply/adjoint use
per-cell **loops** — the asymmetry that makes the 1-D frame three skeletons, not one.

---

## 2. The two shared-frame mechanisms — how they differ

### 2a. `_OctantWalk` (multi-D) — solve+apply sharing by CALLBACK INJECTION ✅ unified

`_interior_walk` (:753) is the shared octant skeleton. It takes **four callbacks** —
`inflow_of`, `shed`, `pure_z`, `interior` — and runs the invariant frame: project the
quadrature octant label to in-plane signs, branch the pure-z degenerate octants, derive the
per-axis in/out domain faces (`_inflow_faces`/`_outflow_faces` :573/:588), read the octant's
inflow, run `interior`, shed each outflow capture. Both public directions route through it.

**Direction-carrying OBJECTS (never a boolean):**
- `_ApplyOperands` (:598) — apply-direction problem data: `probe` (ψ̄), `sig_t`, `str_axes`, `Q_zero`.
- `_SolveOperands` (:619) — solve-direction problem data: `Q` (source), `sig_t`, `str_axes`.
- `_SweepEmit` (:659) — solve OUTPUT mode (angular XOR moment); construction guard (:688) rejects a mixed/empty mode (illegal-state-unrepresentable).
- The "NEVER a boolean `is_solve` flag" tripwire (:738 docstring) is enforced by an **AST tripwire** in `test_one_octant_walk.py:85` (`test_octant_walk_is_kernel_parameterized_not_boolean`): it `ast.parse`s `_OctantWalk`'s source and fails if `is_solve`/`is_apply`/`is_matvec` appear as a real identifier (Name/Attribute/arg — not in docstrings).

**What forks:** ONLY (1) the cell kernel — `cell_kernel_batch` (solve) vs `residual_kernel_batch`
(apply) — and (2) the emit policy (`_SweepEmit` angular/moment vs the matvec's `(L+C)ψ` bulk +
O.4b boundary defect). **What is shared:** the entire octant skeleton.
- `sweep_group` (:798) drives SOLVE — supplies `_SolveOperands` + `_SweepEmit`, boundary I/O via the **live** `boundary_flux` (reads inflow off the trace, sheds outflow back in — the basis for both Jacobi single-group and the G-S inter-group reflect).
- `loss_action` (:858) drives APPLY — accumulates `LpC`, reads inflow from the **GIVEN** `psi.boundary` with **NO `bc.apply`** (the reflective coupling is the sibling `−B`), captures outflow into `streamed`, then computes the **O.4b active-trace boundary residual** (OUTFLOW slots → `streamed − given`; INFLOW slots → identity `given`).

**LANDMINE:** `_interior_walk` has **no transpose**. The adjoint apply is not expressible through
this frame today — it is exactly the missing capability the carve must add.

### 2b. `_OneDimScanWalk` (1-D) — `sweep` and `loss_action` are SEPARATE methods ❌ not unified

CONFIRMED: `sweep` → `_run` (:2816) and `loss_action` → `_apply_walk` (:2250) are **separate
methods**, NOT a shared frame like `_interior_walk`; `loss_action_transpose` (:2602) is a third
standalone body.

- **`_run` (solve)** is a vectorized **Blelloch prefix scan**: `ordinate_scan` over the chain
  order `geom.chain_idx` (:2998), driven by the precomputed two-stratum cache
  (`coll.a_attenuation`/`inverse_denom`/`face_blend_weight`). Slab = exactly **2 scan calls per
  sweep** (one per direction, joint-batched over K ordinates × ng); curvilinear = one scan per
  ordinate per μ-level (the M-M angular thread couples ordinates sequentially).
- **`_apply_walk` (apply)** is a per-cell Python **loop** `for i in cell_indices` (:2408) over
  `dag_walk_cell_indices`, **recomputing** coefficients per cell via `residual_kernel_batch`
  (Cartesian, :2436) / `cell_balance_for_streaming` (curvilinear, :2469).

**Are `_run` and `_apply_walk` true twins?** No — they share **helpers** but not the **walk
skeleton**. Shared: `pole_angular_closure` (the M-M thread + Carlson coupled-pole seed),
`frame_signs_for` (the moment-frame sign involution), and the underlying scheme cell algebra.
Different: scan-over-chain-via-cache vs loop-over-DAG-recomputing. They are structurally
distinct control flows.

**Docstring intent vs reality.** The class docstring (:2153–2155) states the intent:

> "The MATVEC (apply direction) attaches in Phase C as the per-ordinate apply-kernel (the
> α=−1, β=2ψ̄ scan), **mirroring `_OctantWalk`'s cell-kernel injection**."

This was **NOT realized**. `_apply_walk`'s own docstring (:2255) says it was "relocated
**verbatim** off `_MSpatialOperatorSum._compute_LpC`" — a transplanted procedural body, not a
kernel injected into a shared `_run` frame. (Note also: the *intended* α=−1/β=2ψ̄ reflection
scan is what `ScanMarch` 2-D uses (`reflect_scan_coefficients`), but the 1-D `_apply_walk`
Cartesian path does **not** scan — it loops cell-by-cell calling `residual_kernel_batch`.)

**WHY the 1-D walk does NOT use the `_interior_walk` pattern.** The 1-D **solve** is a
parallel-prefix **scan** (a vectorized closed form, structurally a different control primitive
from a per-cell/per-octant DAG walk); `_interior_walk` is a per-cell octant frame. A scan and a
cell-loop cannot share one skeleton without forcing the solve into a loop (losing the cache
vectorization that exists precisely to avoid the per-cell cost — lessons-L16). This is the
documented **scan-vs-wavefront WRONG-FIT** (AGENT.md durable-shape note). Note the 1-D
**apply** and **adjoint** ARE both per-cell loops over `dag_walk_cell_indices`, so *those two*
could share a frame with each other — it is the solve that is the structural odd-one-out.

---

## 3. The DAG layer — "one DAG, reverse traversal"?

**Partly true — and only for 1-D.**

`dag_walk` (augmented_mesh.py:908) / `dag_walk_cell_indices` (:1039) are **1-D ONLY**. They
raise for 2-D Cartesian (`reduced is None` ⇒ "2-D Cartesian wavefront sweeps use anti-diagonal
scheduling, not per-cell visits", :992-998 / :1065). The cell order is `range(nx)` for μ≥0,
`range(nx-1,-1,-1)` for μ<0 (:1099-1102); cylindrical pure-azimuthal degenerate → `range(nx)`
regardless of sign (:1096).

**Forward vs adjoint over the SAME DAG (the hypothesis, confirmed for 1-D):**
- Forward `_apply_walk`: `cell_indices = list(dag_walk_cell_indices(direction_sign, mu_level_idx))` (:2399), then `for i in cell_indices` (:2408).
- Adjoint `loss_action_transpose`: SAME `dag_walk_cell_indices(direction_sign=s, …)` (:2716), then `for i in reversed(cells)` (:2725).

So at the 1-D level the brief's hypothesis holds: same DAG, **forward iterates, adjoint walks
`reversed`** — NOT μ-reversal (the `direction_sign` picks the physical octant DAG; reversing
the within-octant cell order gives the transpose).

**BUT the cell-ordering sources are FOUR, not one** — map them all:
1. **`geom.chain_idx` / `chain_idx_inv`** — the 1-D **solve** (`_run`, :2998). A precomputed permutation array for the vectorized scan.
2. **`dag_walk_cell_indices`** — the 1-D **apply + adjoint** loops (:2399, :2716). An iterator.
   - (1) and (2) are the **same total order** (both `range(nx)`/reversed) materialized two ways — a Pattern-2 duplication of "the 1-D chain order". *(INFERENCE from matching `range(nx)` logic; not cross-referenced in code.)*
3. **`SweepDependencyGraph`** (`loss_representation/sweep_graph.py`) — the **2-D wavefront** (window + full-field), an **anti-diagonal anti-hyperplane DAG** walked by `walk_windowed`/`walk_full`. A *different* DAG abstraction; `dag_walk` explicitly defers to it for d≥2.
4. **Inline row-march** — `ScanMarch` 2-D builds its own ordering: `x_reverse = sx_eff < 0` (:1700/:1811), `y_rows = range(ny) if sy_eff ≥ 0 else range(ny-1,-1,-1)` (:1710/:1825). Marches rows, scans within each. NOT the DAG, NOT `dag_walk_cell_indices`.

**So "one DAG" is really two DAG implementations (1-D linear chain vs multi-D anti-diagonal
graph) plus the scan permutation plus the row-march.** A unified walk must EITHER bridge the
linear-chain iterator and the anti-diagonal graph, OR accept they stay separate and unify
orientation×kernel *within* each dimensionality band. The multi-D wavefront DAG is **not** the
same abstraction as the 1-D scan order.

---

## 4. The historical WHY (the most important part)

The fragmentation is the residue of **five incremental waves landing along different axes**,
with the one unification campaign (Wave 1) scoped to {forward} × {multi-D} and stopping there.

### Wave 1 — the sweep-strategy carve (S0–S6.9, COMPLETE 2026-06-11, #222)
Source: module docstring (:91-118) + `loss_representations.rst` §History (:1832) + the one-walk
theorem (:1514). It replaced a *scattered procedural* dispatch (`transport_sweep` branching on
`reduced is not None`; the matvec branching in **five** operator gates on `not is_1d`; the
oracle reachable only via hand adapters) with the polymorphic `LossRepresentation`.

Phases (each independently bit-identical-gated): **S1** skeleton (protocol, `_DAGWavefront`, 4
thin-wrapper leaves, `supports`/`default_for`/registry). **S2** matvec side — the 5 gates → one
delegating call. **S3** retires hand adapters → FFW becomes the d-generic oracle. **S4** window
→ `frontier_dim = d−1`. **S5.1** `ScanMarch` opt-in (Fork B1). **S6.2** rename
`SweepStrategy→LossRepresentation`, `residual→loss_action`. **S6.3** the walk moved OFF the
operator (`−C` glue collapsed 5×→1×). **S6.4 ONE WALK** (`_OctantWalk._interior_walk`) +
family-owned DAG cache + `diamond.py` pure kernel pair (`sweep.py` dissolved, `WavefrontFlux`
retired). **S6.5 ONE INSTANCE**. **S6.9** Fork-B2 default flip window→ScanMarch.

S6.4 sub-steps (the `_OctantWalk` docstring :746-748): **(a)** routes window + scan-march
**MATVEC** through `_OctantWalk`; **(b)** brings the **SWEEP** frames in (the one-walk spy test
flips xfail→xpass); **(c)** family-owned per-shape DAG cache; **(d)** folds the full-field
oracle; **(e)** storage adapters retired — lifted from `diamond.py` to the walk layer; **(f)**
`sweep.py` dissolved, orchestration relocated into this module. **All (a)–(f) LANDED** (the
diamond.py "S6.4(e): the storage adapters RETIRED" comment :709 and the live one-walk spy test
confirm).

**Scope of Wave 1:** it unified {forward sweep, forward matvec} for {multi-D} into one walk. It
did **not** touch the adjoint (the transpose did not exist yet) and did **not** pull the 1-D
walk into `_interior_walk` (the 1-D solve is a scan — structurally foreign).

### Wave 2 — the LD / coefficient model (#158 / #239 / #240, 2026-06)
**#158** the coefficient model (`residual_kernel_batch`, `affine_scan_coefficients`). **#239**
the lift (`cartesian_scan_coefficients` / `reflect_scan_coefficients` — the *scheme* owns the
diamond `2`). **#240** spatial moments: D5a single-sources the `2g` fold
(`_cartesian_streaming_diagonal`); **D5b-S3** ships the unified all-d **LD moment matvec** (the
trailing `2^d` spatial-moment axis); **Phase 2 Step B** makes σ **explicit** (the caller
single-sources it). This made the **cell kernel** polymorphic (DD/LD) and σ-single-sourced, but
worked within the existing face/walk structure.

### Wave 3 — Wave O: BC extraction + the adjoint (#208, 2026-06)
**O.4a.2 + O.4b** the **BARE sweep**: the reflective coupling `B·ψ.outflow→inflow` moved OUT of
the sweep into the sibling `−B`; the sweep now reads the given inflow trace directly (the
`bc.apply` is gone — discrete_ordinates.rst :6212). Affine-typed operator algebra
(`FluxDisplacement`, typed residual). **O.2b** the **1-D reverse sweep**:
`loss_action_transpose` landed for **1-D only** (the analytic reverse sweep, the angular
second-triangular-factor adjoint `angular_adjoint`, the Carlson coupled-pole seed adjoint). The
**multi-D Cartesian adjoint was EXPLICITLY DEFERRED** — every multi-D `loss_action_transpose`
raises "O.2b lands the 1-D reverse sweep first".

### Wave 4 — #206 the 1-D relocation
**Phase B**: pure relocation of the free 1-D helpers (`_sweep_1d_unified` / `_ensure_geom_cache`
/ `_ensure_coll_cache` / `_run_1d_sweep`) into `_OneDimScanWalk`, bit-identical — a holder frame
*intended* to mirror `_OctantWalk` but where matvec/transpose attach as separate methods.
**Phase C**: the 1-D matvec (`_apply_walk`, relocated **verbatim** off `_compute_LpC`) + the
transpose (`loss_action_transpose`, relocated **verbatim** off `_compute_LpC_transpose`).

### Wave 5 — #257 S8b / #238 the leaf cleanup
`streaming_action` / `streaming_action_transpose` (the σ-free pure-L leaf, single-sourced
through `loss_action` at σ=0). Retired the `M_spatial` / `M_angular_redist` operator-leaf split
(#238) — the fused `loss_action` is the only path.

### Net: WHY is it fragmented? (2–4 sentences, the honest story)
The forward pair `{sweep, matvec}` *was* genuinely unified for multi-D into one
`_OctantWalk._interior_walk` (Wave 1's "one walk", a real, spy-pinned achievement) — but that
campaign's scope was {forward}×{multi-D} and **stopped there**: the adjoint did not yet exist,
and the 1-D solve is a parallel-prefix scan that does not fit a per-cell octant frame. The
adjoint (`loss_action_transpose`) and the 1-D matvec/transpose arrived **later, in separate
waves (O.2b, #206 Phase C), as verbatim relocations of pre-existing procedural bodies**
(`_compute_LpC` / `_compute_LpC_transpose`), each landing as its own standalone method because
there was no shared frame to inject them into (the 1-D solve is a scan; `_interior_walk` has no
transpose). So the fragmentation is **a real-but-partial unification that was never extended
across the adjoint axis or the 1-D scan band.**

### L21 — where ACHIEVED vs CLAIMED-BUT-NOT
The **L21 invariant** ("matvec ≡ sweep — different applications of ONE operator";
loss_representations.rst :166) is:
- **ACHIEVED** for **multi-D forward**: sweep ≡ matvec via `_OctantWalk._interior_walk`, pinned
  by `test_one_octant_walk.py`'s runtime spy (`test_sweep_and_loss_action_hit_one_octant_walk`)
  + the AST tripwire. The **one-instance theorem** (one representation per operator) is pinned by
  `test_apply_and_solve_share_one_representation_instance`.
- **CLAIMED-BUT-NOT** for:
  1. **the 1-D walk** — `_run` (scan) / `_apply_walk` (loop) / `loss_action_transpose` (reverse
     loop) are **three separate skeletons**; the docstring "mirroring `_OctantWalk`'s
     cell-kernel injection" is unrealized intent.
  2. **the adjoint direction entirely** — no transpose routes through `_interior_walk`; the
     multi-D Cartesian adjoint is DEFERRED (raises); only the 1-D adjoint exists, standalone.
  3. **the fourth face** — `sweep_transpose` `(L+C)⁻ᵀb` does not exist at all.

---

## 5. The cell-kernel layer — how many realizations of the per-cell algebra?

`DiamondDifference` (diamond.py:104) declares **THREE capability groups, and there is NO
transpose group** (the explicit class comment :709-729):

1. **Per-cell reference pair** — `update` (:166, solve) / `residual` (:222, apply) — via
   `cell_balance_terms` (the `CellVisit`/`StreamingTerms` form). The slow reference
   (`scheme.update`) path.
2. **Batched kernel pair** — `cell_kernel_batch` (:386, solve) / `residual_kernel_batch`
   (:456, apply) — both route through the **shared `_cartesian_streaming_diagonal` (:325)**, the
   single source of DD's `2g` fold (:444 / :487). The DAG wavefront family + the scan-apply twin.
3. **Scan-family coefficients** — `affine_scan_coefficients` (:501, the 1-D scan; its **own**
   denom — the curvilinear-general `2|μ|A_down + (ΔA/w)c_out + Σ_t·V`, deliberately NOT folded
   into the Cartesian helper, #242 deferred per the long :584-612 rationale) /
   `cartesian_scan_coefficients` (:634, the 2-D row-march — **via** `_cartesian_streaming_diagonal`)
   / `reflect_scan_coefficients` (:692, the apply reflection α=−1, β=2ψ̄).

**The shared cell algebra (the curvilinear / mask path):**
- `cell_balance_for_streaming` (cell_balance.py:120) — `(denom, numer_upstream)`, geometry-blind
  by interface (the closure supplies the angular terms). Consumed by forward `_apply_walk`
  curvilinear (:2469), the forward **adjoint** for ALL geometries (:2731), and
  `DiamondDifference.residual` (:300).
- `cell_balance_terms` (cell_balance.py:262) — the `StreamingTerms` form (used by
  `DiamondDifference.update`/`residual`). A **mild twin** of `cell_balance_for_streaming`: two
  denom-formula sites (one computes the M-M `c_in`/`c_out` inline :313-314, the other receives
  them from the closure), each "single source" within its layer.
- M-M angular closure (`MorelMontryAngularSweep`, pole_angular_closure.py:488): `cell_contribution`
  (:888) returns the per-cell `(denom_term, upstream_numer_term)` and is **shared forward/adjoint**
  (Pattern 2 — the denom is ψ-independent); `angular_adjoint` (:925) is a **separate
  hand-written reverse-mode adjoint** (the reverse M-M recurrence + Carlson sweep), verified
  bit-for-bit against a dense-probe oracle (`derivations/diagnostics/diag_p42_adjoint_oracle.py`).

**Answer to the single-source question:**
- **Forward solve + forward apply DO build on the SAME cell algebra** — Pattern 2 holds for the
  forward pair. Cartesian: `cell_kernel_batch` ≡ `residual_kernel_batch` share `denom`,
  `couplings`, and the outgoing-face reconstruction (only `psi_avg = numer/denom` vs `residual =
  denom·ψ̄ − numer` differ). Curvilinear: both via `cell_balance_for_streaming` + `cell_contribution`.
- **The adjoint reuses the forward COEFFICIENTS but has NO transpose kernel.** `loss_action_transpose`
  recomputes `denom` through the SAME `cell_balance_for_streaming`/`cell_contribution`
  ("every coefficient is ψ-independent… no twin algebra", :2625-2627) and then **hand-transposes
  the arithmetic inline** (the VJP: :2743-2748) plus the separate `angular_adjoint`. So the
  transpose shares *coefficients*, not a *kernel function*.
- **Asymmetry / landmine #1 — two denom functions, agreement by math not by sharing.** Forward-apply
  Cartesian uses `residual_kernel_batch` (the #158/#240 ÷V coefficient model, `g=|μ|/Δ`), but the
  **adjoint uses `cell_balance_for_streaming`** (the older ×V face-area form) for *all* geometries.
  For slab they produce the same denom **two different ways** (the ×V/÷V convention split) — a real
  seam a unified transpose kernel would have to reconcile.
- **Asymmetry / landmine #2 — the 1-D adjoint appears DD/scalar-ONLY (no LD moments).** *(INFERENCE
  from code shapes, not a docstring claim:)* `loss_action_transpose` allocates `psi_bar =
  np.zeros((ng, N, nx))` with **no `moment_tail`**, and carries **no `_reframe`/`frame_signs`** — so
  it does not carry LD's trailing spatial-moment axis, whereas forward `_apply_walk` does (via
  `residual_kernel_batch` + `moment_tail` + `_reframe`). The LD adjoint looks like an unimplemented
  gap. (Curvilinear is DD-only by design — `LinearDiscontinuous` curvilinear is unpublished — so the
  gap is specifically **LD slab adjoint**.) **Verify before the carve.**
- LD's own multi-D kernel: `LinearDiscontinuous.cell_kernel_batch`/`residual_kernel_batch` are
  d=1-only (raise on d=2; loss_representations.rst :1495). 2-D LD ⇒ wavefront only; `ScanMarch`
  refuses it (the `transverse_coupling_is_facewise` gate).

---

## 6. The landmines for the unification carve

### A. Face → canary map (the bit-identity pins each face leans on)
| Face | Canary(ies) |
|------|-------------|
| `sweep` 1-D solve | `tests/sn/sweep/core/test_sweep_cache.py` (cache dual-view, **rtol=1e-13**) + slab snapshots rtol=1e-12 |
| `sweep` multi-D solve | `window ≡ full` oracle (`test_sweep_graph_window_equivalence`); `tests/sn/solve/test_scan_march_end_to_end.py` (Mode-9 FP-invariance, forced-window gates) |
| `loss_action` matvec | `tests/sn/operators/test_invertible_operator.py` (`test_invertible_apply_is_M_of_C_sigma_bit_identical`); the FFW oracle (G2.c); `tests/sn/sweep/core/test_cell_kernel_batch.py` |
| `loss_action_transpose` adjoint | `tests/sn/operators/test_g_adjoint_reciprocity.py` — `⟨Aψ,φ⟩_G == ⟨ψ,Aᵀφ⟩_G` for slab/sphere/cyl + 2g (`@pytest.mark.foundation`), the **L11 wrong-trace-metric** negative control, and the dense-probe oracle `diag_p42_adjoint_oracle.py`; `tests/sn/sweep/core/test_phase_c_gates.py` |
| one-walk structure | `tests/sn/operators/test_one_octant_walk.py` — runtime `_interior_walk` spy + the AST tripwire (`is_solve`/`is_apply`/`is_matvec`) |

### B. Deferred faces the unified abstraction must cleanly express as "raises here"
- **multi-D Cartesian adjoint**: `_DAGWavefront.loss_action_transpose` (:1091) and
  `ScanMarch.loss_action_transpose` 2-D (:1858) — both `NotImplementedError`. (`_OneDimScanWalk`
  itself also guards multi-D Cartesian out of `_apply_walk`/`loss_action_transpose` at :2316/:2647.)
- **2-D LD**: `ScanMarch` refuses (supports gate); the wavefront LD cell kernel is d=1-only.
- **`sweep_transpose` `(L+C)⁻ᵀb`**: does not exist — no method, no operator door. The unified
  walk that can run a reverse traversal with the solve kernel would, in principle, *be* this
  fourth face — but nothing consumes it today (no Krylov/adjoint-solve caller wired), so adding
  it is genuinely new surface, not just a re-layering.

### C. "Orientation" is MORE than `reversed()` (the curvilinear seams)
For **slab** (Cartesian + `IdentityAngularClosure`) the adjoint **IS** ≈ `reversed(cells)` +
the boundary in↔out swap — the clean 2×2 hypothesis is closest to reality here. Curvilinear
adds two extras the unified "orientation" parameter must carry:
1. **Boundary in↔out swap**: forward reads inflow / writes outflow; adjoint reads OUTFLOW
   cotangents / writes INFLOW cotangents (:2677-2695, :2772-2775).
2. **Carlson coupled-pole MIRROR seed** (curvilinear only): forward runs the −1 sweep **first**,
   and its pole outflow read at the *mirror* ordinate `reflection_index("x")` seeds the +1 sweep
   (:2492/:2505/:2512). The adjoint **reverses this cross-direction data dependency**: the +1
   seed cotangent routes into the −1 reversal's outflow cotangent at `mirror[gd]` (:2758). This
   is not `reversed(cells)` — it is a reversed *cross-sweep* coupling within a μ-level (pinned by
   the μ-level-preservation invariant, discrete_ordinates.rst :8689).
3. **Angular SECOND triangular factor** (curvilinear M-M recurrence in the ORDINATE index): the
   forward sweep is lower-triangular in BOTH cell-visit order AND the ordinate index (two
   triangular factors); the forward builds the angular factor via `precompute_psi_state` +
   `cell_contribution`, the adjoint reverses it via the separate `angular_adjoint` (:2765), zero
   for slab.

### D. The scan-vs-loop EXECUTION gap (the deepest landmine)
The **solve** is a vectorized prefix scan (2 `ordinate_scan` calls/sweep for slab, driven by the
precomputed cache); the **apply/adjoint** are per-cell Python loops (recomputing coefficients).
Unifying solve+apply into one walk forces a choice:
- force the solve into a per-cell loop ⇒ **loses the scan vectorization** the cache exists to
  provide (a measured perf regression; lessons-L16); OR
- force the apply into a scan ⇒ the α=−1/β=2ψ̄ reflection-scan exists (`reflect_scan_coefficients`)
  and **`ScanMarch` 2-D already uses it for the x-row**, but the **1-D `_apply_walk` does not**
  (it loops). Even within "apply" the vectorization strategy is inconsistent between 1-D and 2-D.

### E. The DAG is two abstractions (re §3)
`dag_walk_cell_indices` (1-D linear chain) vs `SweepDependencyGraph` (multi-D anti-diagonal) vs
the `geom.chain_idx` scan permutation vs the `ScanMarch` inline row-march. A literal "ONE
per-ordinate cell DAG" across all dims does not exist today.

---

## What makes the carve HARDER than the clean orientation×kernel 2×2

The clean hypothesis is **{forward, reverse} × {solve, apply}**. Reality adds:
1. **A third, non-free axis — execution strategy {scan / cell-loop / wavefront-graph}.** The
   solve *needs* the scan (perf); the multi-D path *needs* the anti-diagonal graph; these are not
   interchangeable with the per-cell loop the apply/adjoint use.
2. **"Orientation" is not a clean involution in curvilinear.** The Carlson mirror seed couples the
   ±1 sweeps cross-wise, and the angular factor is a *second* triangular structure reversing
   independently. Only slab is the clean `reversed()` + swap.
3. **The grid is asymmetric in coverage.** `sweep` (all) / `loss_action` (all) /
   `loss_action_transpose` (1-D only) / `sweep_transpose` (nonexistent). The unified walk must
   express a partial/deferred grid with the existing `NotImplementedError` deferrals intact.
4. **The LD asymmetry** (forward-apply carries LD moments; the 1-D adjoint appears not to —
   verify) and the **two-denom-function seam** (×V `cell_balance_for_streaming` vs ÷V
   `residual_kernel_batch`).

**The cleanest sub-target** (where the 2×2 is closest to true): **slab/Cartesian**, where the
adjoint already *is* `reversed(cells)` + boundary swap over the SAME `dag_walk_cell_indices`,
and the forward solve/apply/adjoint already share `_cartesian_streaming_diagonal`'s coefficients.
The 1-D apply and adjoint (both per-cell loops over `dag_walk_cell_indices`) are the natural
**first** unification (orientation = forward-iterate vs `reversed`, kernel = apply vs the
hand-transposed apply) — the solve-scan and the multi-D adjoint are the genuinely hard parts.

---

### Provenance / confidence notes
- **STATED (read directly):** the one-walk theorem scope, all deferrals/raises, the "verbatim
  relocation" provenance, the three DiamondDifference capability groups + no-transpose-group,
  Pattern-2 sharing of the forward pair, the Carlson mirror seed + `angular_adjoint`, the DAG
  1-D-only restriction.
- **INFERENCE (flagged inline):** (a) the 1-D adjoint being DD/scalar-only with no LD-moment
  support — read from buffer shapes in `loss_action_transpose`, not a docstring; (b)
  `geom.chain_idx` ≡ `dag_walk_cell_indices` total order — from matching `range(nx)`/reversed
  logic, not a code cross-reference. Both should be confirmed before the carve commits to them.
