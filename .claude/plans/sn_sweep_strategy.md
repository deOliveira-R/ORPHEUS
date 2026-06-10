# SN Sweep-Strategy Carve (C3.4 + C3.5) — polymorphic sweep/matvec dispatch

## ⭐ STATUS / RECOVERY (read first)

Part of the N-D layout campaign. **Authoritative index = `.claude/plans/nd_layout_foundation.md`**;
this file is the self-contained design for the sweep-strategy carve, which **re-scopes** the
original C3.4 ("wire a `_wavefront_1d_sweep` adapter") + C3.5 ("orchestration d-generic") around a
first-class `SweepStrategy` abstraction.

- **Designed 2026-06-10** in a multi-turn design conversation; every decision below is LOCKED
  (see §"Decisions locked"). **NOT yet implemented.**
- **Depends on C3.0–C3.3 (all DONE):** the wavefront spine is already dimension-generic — the
  per-octant DAG `SweepDependencyGraph.from_cartesian(shape)` (C3.1), the diamond cell kernel
  `cell_kernel_batch` (C3.2a), the full-field walk `graph.apply`/`graph.residual` (C3.2b), and the
  geometry octant build + `is_1d`/5-gate dispatch (C3.3). The spine is B7-verified at d=1 **and**
  d=3 (`test_sweep_graph_nd_admission`).
- **First implementation step is a `test-architect` dispatch** — this carve routes the *matvec*
  from one shape to another across a subsystem boundary, which is the proactive trigger in
  `subagent-handoff-protocol` ("operator-algebra carve crossing subsystem boundaries"). Do NOT
  start S2 before the verification plan lands.
- **Branch:** `worktree-sn-nd-layout` (off `main`); 33 commits ahead of origin, NOT pushed.

---

## THE PROBLEM THIS SOLVES

The 1-D-vs-multi-D sweep dispatch is scattered and procedural, and it is the *same decision*
spelled three different ways:

1. **`transport_sweep`** (`orpheus/sn/sweep.py:104`) branches on dimensionality:
   1-D → `_sweep_1d_unified` (the Blelloch `ordinate_scan`), 2-D → `_sweep_2d_wavefront`
   (the `_MovingFrontier` window).
2. **The matvec** branches in **5 operator gates** (C3.3: `not sn_mesh.is_1d` in
   `_compute_LpC`, `_compute_decomposition`, `_compute_LpC_transpose`, `apply`,
   `apply_transpose`).
3. **The verification oracle** (the full-field spine `graph.apply`) is reached only through
   **hand-built TEST adapters** — `_wavefront_1d_sweep` and `_cumprod_1d_sweep` in
   `tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py` — because the spine is not a
   *selectable* path. There is no way to ask "which sweep methods are valid for this mesh?".

Adding a method, a dimensionality, or a frontend means touching all three. That is an
enum-style branch repeated at every call site — **cyclomatic complexity, not abstraction**
(`coding-elegance` anti-pattern: stringly-typed / special-case dispatch). An enum parameter
threaded into `transport_sweep` would only make it worse (a second branch axis).

---

## GOVERNING PRINCIPLE (the rule for every strategy)

> **Construct each strategy as general as its algorithm naturally allows. Select narrow.
> Specialize the implementation only on a *measured internal* performance cost.**

Three separable layers, never conflated:

- **Construct general (capability).** If the algorithm is *naturally* d-general, the CODE handles
  any d. A strategy that is d-specific in code *only because its algorithm is intrinsically
  d-specific* is legitimate (a prefix scan needs a total order → a chain → 1-D; there is no
  "2-D chain scan"). A strategy that is d-specific *because we wrote a narrow crutch* is a smell —
  it means the general structure did not surface.
- **Select narrow (policy).** Whether we OFFER / recommend / default a strategy at a given
  `(geometry, ndim)` is a SEPARATE layer (`supports` / `default_for`), independent of the code's
  capability. "Don't pick the window at d=1, pick cumprod" is a recommendation, *not* a reason to
  leave the window's code unable to express d=1.
- **Specialize on measured cost only.** The sole justification to restrict an implementation's
  d-range is a *measured* hot-path regression where the general construction makes the *kept*
  cases slower than a narrower construction would. Even then, keep the general path as the pinned
  fallback/oracle (the `feedback_aggressive_retirement` "fuller-view oracle" exception).

This principle is the acceptance lens for every phase below.

---

## THE ARCHITECTURE

### The strategy protocol (the algebra)

```python
class SweepStrategy(Protocol):
    def sweep(self, Q, sig_t, boundary, *, initial_guess=None) -> tuple[...]:  ...  # (L+C)⁻¹ q
    def residual(self, psi, ...) -> ...: ...                                        # (L+C) ψ   (matvec twin)
    @classmethod
    def supports(cls, mesh) -> "Compatibility": ...                                 # selection layer
```

Each strategy carries **both** the forward sweep AND the matvec twin, so the operator's `apply`
routes through `strategy.residual(...)` and the 5 C3.3 gates collapse.

### The hierarchy (what is shared vs not)

```
SweepStrategy (Protocol: sweep, residual, supports)
├── _DAGWavefront            ── reads mesh.sweep_graphs (the anti-hyperplane DAG) + the DD kernel
│   ├── FullFieldWavefront     buffer = full field     · Cartesian, any d · the ORACLE
│   └── MovingFrontierWindow   buffer = rolling (d−1)-frontier · Cartesian, d≥2 (built general) · prod opt
└── CumprodScan             ── reads the two-stratum scan cache · 1-D, any geometry · prod opt
```

**Substrate lives on the mesh, not the strategy** (strategies are lightweight consumers):
- `sweep_graphs` — the per-octant anti-hyperplane DAG (`from_cartesian(spatial_shape)`), built only
  for Cartesian meshes (C3.3); curvilinear sets it `None`.
- the two-stratum scan cache (`geom`+`coll` → `ordinate_scan`) — the chain-recurrence substrate.

**Is it the same DAG?** Yes for two of three: `FullFieldWavefront` and `MovingFrontierWindow`
consume the **same** `sweep_graphs` — they are two *buffer policies* over one anti-hyperplane walk
(full field vs rolling `(d−1)`-frontier), already pinned bit-identical by the C3.2b
`window ≡ full` oracle. `CumprodScan` builds no DAG — 1-D is a chain, the Blelloch closed form
needs no graph. The DD physics (closure, attenuation) is shared *conceptually* across all three but
expressed two ways (scan recurrence vs explicit kernel), pinned equivalent by the equivalence tests
— the accepted dual-view.

### Why these three (per the governing principle)

- **`CumprodScan`** — intrinsically 1-D (chain prefix scan needs a total order). Geometry-blind:
  slab, sphere, cylinder share one body via `_sweep_1d_unified` → `ordinate_scan` (the angular
  redistribution folds into the scan's affine source `b`, ordinates processed in angular order).
  Legitimately 1-D *by the algorithm's nature*, not by our hand.
- **`FullFieldWavefront`** — naturally d-general; ALREADY d∈{1,2,3} (the DAG + kernel + walk are
  d-generic, B7-verified d=1/d=3). The slow, readable spine; the verification oracle.
- **`MovingFrontierWindow`** — naturally d-general (the moving frontier is the `(d−1)`-dim rolling
  slab: a point at d=1, a line at d=2, a surface at d=3). **Must be constructed on
  `frontier_dim = d−1`**, capable of any d (S4). The 2-diagonal is the `frontier_dim == 1`
  instance, not the thing itself.

---

## SELECTION / COMPATIBILITY MODEL (the frontend-checkable layer)

Applicability is a **declared, queryable capability** — "make illegal states unrepresentable"
applied to method selection. The compatibility signal is the genuine criterion (the coordinate
system), NOT the `sweep_graphs is None` substrate proxy:

```python
# mesh.is_cartesian  ::=  curvature is None   (new one-line property, parallel to is_1d)

class CumprodScan:
    @classmethod
    def supports(cls, mesh): return Compatibility(mesh.is_1d, "requires a 1-D mesh")

class FullFieldWavefront:
    @classmethod
    def supports(cls, mesh): return Compatibility(mesh.is_cartesian, "requires Cartesian geometry")

class MovingFrontierWindow:
    @classmethod
    def supports(cls, mesh):
        return Compatibility(mesh.is_cartesian and mesh.ndim >= 2,
                             "requires Cartesian geometry, d ≥ 2")
```

`Compatibility(ok: bool, reason: str)` — the `reason` lets a teaching frontend gray-out a method
*and explain why* ("Moving-frontier window — requires Cartesian geometry, d ≥ 2"), which is
pedagogically load-bearing (ORPHEUS teaches reactor physics).

**Three consumers, ONE predicate (single source of truth):**
1. **Frontend** — `[S for S in SWEEP_STRATEGIES if S.supports(mesh).ok]`. A cylinder
   (non-Cartesian) → only `CumprodScan` → the dropdown shows only Blelloch.
2. **Factory default** — `SweepStrategy.default_for(mesh)` picks the best *available* production
   optimization, **falling back to the spine** when no optimization exists:

   | mesh | applicable | `default_for` |
   |---|---|---|
   | Cart-1D | `{FullField(oracle), CumprodScan}` | `CumprodScan` |
   | Cart-2D | `{FullField(oracle), MovingFrontierWindow}` | `MovingFrontierWindow` |
   | Cart-3D | `{FullField}` (window not built yet) | `FullField` ← never stuck |
   | Cyl/Sph-1D | `{CumprodScan}` | `CumprodScan` |

   The day a d=3 window lands (S4 + a 3-D mesh), Cart-3D's default flips from `FullField` to the
   window automatically — one predicate, no caller change.
3. **Construction guard** — `Strategy(mesh)` raises `IncompatibleStrategy` if `not supports.ok`, so
   even a bypassed UI cannot build an illegal pairing.

That `supports` predicate **is** the `is_1d`/`curvature` dispatch scattered across `transport_sweep`
+ the 5 operator gates today — now declared once per strategy. The whole point: "add 3-D window
support" becomes "extend one strategy, widen one predicate," not a hunt through every call site.

---

## MATVEC UNIFICATION (the C3.3-gate collapse)

The 5 C3.3 gates each spell `not sn_mesh.is_1d`. With `strategy.residual`, they collapse to
"ask the mesh for its strategy and delegate":

- `StreamingOperator.apply` (the live 2-D dispatch, `operator.py:~1458`) → `strategy.residual(...)`.
- The 4 raise-gates (`_compute_LpC`, `_compute_decomposition`, `_compute_LpC_transpose`,
  `apply_transpose`) → the strategy either implements the path or its absence is the
  `IncompatibleStrategy`/not-applicable signal (no hand-written "multi-D not wired" raise — the
  strategy's `supports` carries that).
- The operator holds `strategy = SweepStrategy.default_for(sn_mesh)` (selected once at
  construction) and the hot path is branchless `strategy.residual(...)`.

This is the operator-algebra carve crossing a subsystem boundary → **`test-architect` proactive
dispatch is mandatory before S2.**

---

## VERIFICATION STRATEGY

- **Bit-identical default gate (S1–S2).** `default_for(mesh)` must pick the SAME path the old
  dispatch did, so `transport_sweep` and the operator matvec are byte-for-byte unchanged on every
  existing mesh. The strategies in S1/S2 are *thin wrappers over the existing code* — pure
  dispatch refactor, zero algorithm change. Gate: the C3.2b/C3.3 bit-identity anchors
  (`test_affine_carve_bit_identity`, the DD regression snapshot, A2D-1 source hash) stay green.
- **Solve-vs-solve equivalence (S3).** Replace the hand-built adapters with strategy-vs-strategy
  through the real API: `CumprodScan(mesh).sweep(...)` vs `FullFieldWavefront(mesh).sweep(...)` at
  d=1; `MovingFrontierWindow(mesh)` vs `FullFieldWavefront(mesh)` at d=2; folding in
  `test_2d_full_field_oracle`. ONE pattern, parametrized by mesh ndim.
- **vv-principles Mode 9 (splitting/acceleration FP-invariance).** The equivalence tests run on an
  **anisotropic / heterogeneous** config (NOT the degenerate isotropic-reflective box), asserting
  the converged value is strategy-invariant to solver tolerance — a swapped optimization MUST NOT
  move ψ*. The analytical anchor (`k_inf = 1.875`, `test_cumprod_path_hits_analytical_kinf`) pins
  the converged value so the spine inherits it.
- **Principled-equivalence, not bit-identity, across strategies.** Cumprod vs spine differ at
  FP-association ULP (`assert_array_almost_equal_nulp`) — the existing `_NULP_BOUND` gate.
- **Synthetic d=3 admission (S4).** The general `MovingFrontierWindow` is verified `window ≡ full`
  on a synthetic 3-D shape (B7-style, no 3-D quadrature) to PROVE the `frontier_dim` construction is
  genuinely general — correctness, separate from the perf-quality question.

---

## PHASES (each independently bit-identical-gated + committable)

**S0 — `test-architect` dispatch (PROACTIVE, mandatory).** Verification plan for the carve:
which existing tests pin the legacy dispatch, the solve-vs-solve equivalence specs, the
bit-identity default gate, the vv-Mode-9 config, the compatibility/registry/guard tests. Output
shapes S1–S5.

**S1 — Strategy skeleton + 3 strategies as thin wrappers, SWEEP side.** `SweepStrategy` protocol +
`_DAGWavefront` base + `FullFieldWavefront`/`MovingFrontierWindow`/`CumprodScan` wrapping the
EXISTING sweep code (no behavior change) + `is_cartesian` property + `supports`/`default_for` +
`SWEEP_STRATEGIES` registry. Rewire `transport_sweep` → `default_for(mesh).sweep(...)`.
**Bit-identical** (each strategy wraps the path the old branch chose).

**S2 — Matvec side: `strategy.residual` + collapse the 5 gates.** Each strategy gains `.residual`
(wrapping `residual_windowed` / `graph.residual` / the 1-D matvec). Route `StreamingOperator.apply`
+ the 4 raise-gates through the strategy. **Bit-identical** (wraps existing matvec paths). Operator
holds `default_for(sn_mesh)`. *(Operator-algebra carve — S0 test-architect plan governs.)*

**S3 — Solve-vs-solve equivalence + retire adapters.** Rewrite `test_wavefront_cumprod_equivalence`
+ fold `test_2d_full_field_oracle` as strategy-vs-strategy, parametrized by ndim, on an
anisotropic/heterogeneous config (Mode 9). **Retire** `_wavefront_1d_sweep`, `_cumprod_1d_sweep`,
and the per-d oracle drivers `_sweep_2d_full_field` (`feedback_aggressive_retirement` +
`feedback_retirement_means_test_migration`). Measure the cumprod-vs-spine speedup through the
selector (the original C3.4 perf ask).

**S4 — `MovingFrontierWindow` `frontier_dim = d−1` generalization (the principle, embodied).**
Refactor `_MovingFrontier` / `_window_metadata` from hardcoded 2-diagonal to the general
`(d−1)`-frontier rolling slab. Verify d=1 (the trivial `frontier_dim == 0` base — free, no cost to
d≥2: a 0-dim frontier is a clean degenerate base, no hot-loop branch) + d=2 (bit-identical) +
**synthetic d=3 admission** (`window ≡ full` on a synthetic shape). **FLAG, do not assume:** the
d=2 speedup is a *zero-copy contiguity* win (the frontier *line* stored compactly is basic-slice
addressable; the full-field grid-diagonal needs fancy-indexed copies). Whether the general
`frontier_dim` layout preserves that contiguity at the d=3 *surface* is the ONE place the governing
principle's "measured cost" exception may bite — settled by profiling the d=3 frontier, not assumed.
If it loses the speedup, the fix is a measured d=2 contiguous *fast-path kept alongside* the general
frontier (pinned equivalent), never a d=1 exclusion.

**S5 — Frontend compatibility + cleanup + docs.** `Compatibility(ok, reason)` finalized for the
frontend. Retire the C3.2b elegance CONCERN — the `str_axes` hand-listed axis tuple at the
orchestrators → the strategy holds an `axes`-keyed `sn_mesh.streaming(a)` map (ONE axis-order
source). Retire the `OctantLabel.sign_x/y` / `streams_in_2d` 2-D shims. **Sphinx docs** for the
strategy architecture (dispatch the `archivist`): the protocol, the hierarchy, the compatibility
grid, the governing principle, the spine-as-always-available-capability story.

*Mapping to the original plan:* S1+S2 = the "C3.4/C3.5 strategy + matvec together" the user chose;
S3–S5 subsume the original C3.4 (cumprod oracle + speedup) and C3.5 (orchestration d-generic,
str_axes fix, shim retirement). The original C3.6 (end-to-end 3-D — needs a 3-D mesh/quadrature)
remains downstream.

---

## DECISIONS LOCKED (2026-06-10 design conversation)

1. **Two selectable wavefront strategies** (`FullFieldWavefront` + `MovingFrontierWindow`) sharing
   a `_DAGWavefront` base — NOT one strategy with a buffer-policy flag. Rationale: the frontend
   should let a user pick the oracle explicitly for verification, and they have genuinely different
   memory profiles (the 3× peak-memory split C3.2b measured).
2. **The wavefront oracle is Cartesian-only.** Curvilinear is not an anti-hyperplane lattice; its
   verification reference is the per-cell `cell_update` march the scan is already pinned against.
   No apples-to-oranges oracle. Curvilinear meshes keep `sweep_graphs = None`.
3. **Selection keys on geometry first, via `is_cartesian`** (the genuine criterion), not the
   `sweep_graphs is None` substrate proxy. Honestly a 2-axis `(coord × ndim)` compatibility; the
   DAG family is gated on Cartesian, `CumprodScan` on the orthogonal `ndim == 1`.
4. **Strategies are constructed as general as their algorithm naturally allows** (the governing
   principle). `FullFieldWavefront` already d∈{1,2,3}; `MovingFrontierWindow` built on
   `frontier_dim = d−1` (any d); `CumprodScan` 1-D by the algorithm's nature, not by our hand.
5. **Selection ≠ implementation generality.** "Don't recommend the window at d=1" lives in
   `supports`/`default_for`; it does NOT justify a d=2-hardcoded window. The only implementation
   d-restriction allowed is a *measured* internal performance cost.
6. **The general path is never absent.** `FullFieldWavefront` (the spine) is the applicable
   fallback for any Cartesian d, so `default_for` is never stuck (Cart-3D defaults to the spine
   until a window is built). The code is never "incapable" at any Cartesian dimensionality.

---

## RISKS / OPEN

- **d=3 window contiguity (S4).** The one genuine "measured performance price" candidate — the
  zero-copy speedup may be d=2-line-contiguity-specific. Decided by profiling, not assumed. Memory
  win generalizes for free; speed win is the open question.
- **Bit-identity discipline (S1–S2).** The whole value of S1/S2 is that they are pure dispatch
  refactors. If any strategy's wrapper diverges byte-for-byte from the path it replaces, the carve
  has a hidden behavior change — the bit-identity anchors are the tripwire.
- **`residual` surface area (S2).** The matvec twin must cover forward + adjoint
  (`apply_transpose`) for each strategy. The curvilinear adjoint + the 2-D adjoint are partly
  deferred today (the C3.3 raise-gates) — the strategy must preserve those deferral boundaries
  (an unsupported adjoint is `IncompatibleStrategy`/not-applicable, not a silent wrong answer).
- **Frontend is hypothetical.** No frontend exists yet; `supports`/`Compatibility` are built for
  the factory + tests now, frontend-ready by construction. Do not build UI here.

---

## VERIFICATION GATES / ACCEPTANCE

- S1–S2: every existing sweep + operator + eigenvalue + solve + curvilinear suite green under
  `python -O`; the bit-identity anchors (`test_affine_carve_bit_identity`, DD regression snapshot,
  A2D-1 source hash) UNCHANGED; elegance-enforcer PASS (no enum/branch dispatch survives in a hot
  path; selection is single-sourced in `supports`/`default_for`).
- S3: solve-vs-solve equivalence green at the `_NULP_BOUND` on an anisotropic/heterogeneous config;
  `k_inf = 1.875` anchor intact; zero remaining `_wavefront_1d_sweep`/`_cumprod_1d_sweep`/
  `_sweep_2d_full_field` references (retirement audit).
- S4: `window ≡ full` bit-identical at d=2 (unchanged); d=1 + synthetic d=3 `window ≡ full` green;
  the speedup measurement recorded (d=2 confirmed; d=3 flagged pending or measured).
- S5: Sphinx builds clean; the `str_axes` + `OctantLabel` shims retired (retirement audit);
  no orchestrator hand-lists a positional per-axis tuple.
