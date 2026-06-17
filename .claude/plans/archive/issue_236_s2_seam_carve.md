# #236 §2 — unified CellUpdate-seam refactor (every strategy consumes the scheme)

**Status (updated 2026-06-14, session 3): capability + cache delegation DONE + verified. PIVOTED to #206.**
The session-3 dependency audit revealed the 1-D helpers (`_run_1d_sweep`/`_sweep_1d_unified` for the
sweep + the operator's `_compute_LpC`/`_compute_decomposition` for the matvec) are a **SHARED 1-D scan
core** consumed by BOTH `CumprodScan` AND `ScanMarch` (its `s_y=0` 1-D degeneration; `ScanMarch.supports`
= `is_1d or 2-D-cartesian`) — NOT CumprodScan-private legacy. So "full unification (symmetry with the
wavefront family)" is NOT "fold the helpers into CumprodScan" (that would strand ScanMarch's 1-D path);
it is a **shared-1-D-scan-frame extraction + moving the matvec off the operator** (the unfinished 1-D
half of S6.3b `3a79ab3`, entangled with `_compute_decomposition` → transient orchestrator / `apply` /
`M_angular_redist`) = **#206** ("Unify matvec/sweep WDD recurrence machinery"). The user (session 3)
chose to pivot to the full #206 carve now. **⟹ The carve plan is now `.claude/plans/issue_206_carve.md`
(READ FIRST). This file's closure-routing (old steps 3–5) becomes #206 Phase A; its capability design
(steps 1–2) is DONE.**
Created 2026-06-13 session 2. Parent: `.claude/plans/sn_space_angle_discretization_plan.md`.
Supersedes the earlier "curvilinear-only" framing (user steer: ALL strategies must consume the
selectable scheme).

## ⭐⭐ DONE so far (session 3, on branch `refactor/sn-cellupdate-seam-slab`, uncommitted)

- **Step 1 ✅ — the `CellUpdate` seam capability.** `is_affine_scannable` trait on the `CellUpdate`
  Protocol + `CellUpdateBase` (`=False` default) + `DiamondDifference` (`=True`); three DD methods —
  `affine_scan_coefficients` (σ_t-epoch `(a_attenuation, inverse_denom)`, math lifted verbatim from the
  old `CollisionCache.from_geometry` op-order — the TRAP satisfied), `cell_average_from_faces` (diamond
  `½(in+out)`), `outgoing_face_from_average` (diamond `2ψ̄−in`). NotImplementedError defaults on Base
  (mirroring `cell_kernel_batch`). Two test-mock conformers got `is_affine_scannable=False`.
  ⚠ NAMING: I used domain names `cell_average_from_faces`/`outgoing_face_from_average` (NOT the plan's
  `scan_closure`/`transverse_close`) — name by what-it-is.
- **Step 2 ✅ — cache delegates to the scheme.** `CollisionCache.from_geometry(geom, sig_t, cell_update)`
  now calls `cell_update.affine_scan_coefficients(...)`; `cumprod_a` + the sig_t chain-reorder stay in
  the cache. Threaded `cell_update` through all 4 production callers (solver ×2, operator
  `materialize_inverse_cache`, `_ensure_coll_cache`) + 2 test files. **Producer-side fix (L18):** retyped
  `SNMesh.cell_update` (Protocol `CellUpdate`) → `CellUpdateBase` (the optional capabilities live on Base;
  every scheme IS a Base subclass) — cleared 4 baseline #226 errors at `_CellSolve`/`_CellResidual`.
- **Verified:** bit-identity probe (slab/sphere/cyl, multi-group) + `test_dd_regression` (full solves
  through the cache) + `test_sweep_cache` + sweep/core 437 + cartesian_2d + curvilinear streaming-equil
  ALL GREEN. Pre-existing reds (NOT mine, confirmed via stash on clean HEAD): `test_sphere_{1g,2g}_apply_bit_identical`
  (stale post-ERR-058 snapshots) + #212 `test_heterogeneous_absolute_keff` hang — both curvilinear/eigenvalue, orthogonal.

## ⭐ START HERE (cold-start) — SUPERSEDED by #206; read `.claude/plans/issue_206_carve.md` first

**Steps 1–2 are DONE (see the DONE block above).** The remaining steps 3–7 below are the
**closure-routing**, now reframed as **#206 Phase A** (route the SHARED 1-D scan core's sweep +
matvec closures through `cell_update`). The detail below is still accurate as #206-Phase-A reference;
the surrounding #206 carve (shared-frame extraction + matvec off the operator) is the new scope. The
original "BUILD PR1 slab" framing is retained below for the closure-routing mechanics only.

- **Branch:** already on `refactor/sn-cellupdate-seam-slab` (off `main` @ `eab05ab`); **EMPTY — no
  commits, no production code written yet.** If not on it: `git checkout refactor/sn-cellupdate-seam-slab`
  (it exists) — else recreate off `main`.
- **Done in this thrust:** Tier 0 #234 (merged `eab05ab`, on `main`); explorer seam-map + unified-
  refactor study; test-architect verification spec; ALL architecture decisions (this file). `main` green.
- **Ground from:** this whole file + parent `.claude/plans/sn_space_angle_discretization_plan.md`
  (top Status block) + cluster memory `project-curvilinear-sn-cluster` (auto-loaded). The empirical
  probe is `$CLAUDE_JOB_DIR/tmp/probe_cellupdate_seam.py` (may be gone post-compaction — it's a
  promotion candidate → `test_seam_honored.py`; re-derive from §"Verification" if needed).
- **NEVER `git add`:** `.claude/plans/r1_phase_a_dim_agnostic_ultraplan.md`,
  `derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py`, `scratch/literature/`. Stage
  explicit paths only. Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
  Host env → `.venv/bin/python`; canonical test `python -O -m pytest`; bg tmp = `$CLAUDE_JOB_DIR/tmp`.

### PR1 — ordered steps
1. **Add the capability** to `CellUpdate` Protocol + `CellUpdateBase` (NotImplementedError defaults,
   mirroring `cell_kernel_batch` at `cell_update.py:573`): `is_affine_scannable: bool` (ClassVar),
   `affine_scan_coefficients(...)`, `scan_closure(...)`, `transverse_close(...)`. Implement all on
   `DiamondDifference` (`is_affine_scannable=True`). **← SURFACE THIS FOR USER REVIEW before step 2.**
2. **Rewire** `CollisionCache.from_geometry` (`sweep_cache.py:393`) → call
   `DD.affine_scan_coefficients` for `(a_attenuation, inverse_denom)`; keep `cumprod_a=cumprod(a)` in
   the cache. `test_cache_populator_matches_cell_balance_terms` stays green (now a true contract).
3. **Route the slab closure** `_run_1d_sweep:2004` (`0.5*(in+out)`) → `cell_update.scan_closure`.
4. **Route the slab matvec** `operator.py _compute_LpC:469` (`2.0*psi_cell - psi_face_in`) → the
   cell-update apply form (the closure's apply direction / `residual`).
5. **Gate the scan schedules:** add `and mesh.cell_update.is_affine_scannable` to `CumprodScan.supports`
   (`:701`) + `ScanMarch.supports` (`:1213`). (Slab+non-affine then falls through to `FullFieldWavefront` ✓.)
6. **Tests** (test-architect spec): `tests/sn/regression/test_sweep_seam_carve.py` (single-sweep
   `--capture-baseline`, `kind="direct"`, `array_equal` rtol=0 — A1/A3); `tests/sn/sweep/core/test_seam_honored.py`
   (promote the probe; +/− controls; `np.testing`/`pytest.fail` for `-O` — B1/B2/B3).
7. **qa + elegance-enforcer**; full `-O` slab + regression + sweep-core suites green; `python -m
   tests._harness.audit` exit 0; Sphinx clean; **ff-only merge** to `main`; a `chore(claude)` commit
   captures this plan + the on-disk agent memos with the merge.

## The finding (study-confirmed, structural + runtime)

Of the four SN sweep strategies, the two **production defaults** bypass the `CellUpdate` scheme;
the window + oracle route through it:

| Strategy | closure mechanism | routes through scheme? |
|---|---|---|
| `CumprodScan` (1-D) | inlined `ordinate_scan` + `½(in+out)` (`:2004`/`:2167`); coeffs from the L15 cache | ✗ |
| `ScanMarch` (2-D) | inlined per-row `(α,β)` + `½(in+out)`/`2ψ̄−in_y` (`:1301-1310`,`scan.py:331-332`,`:1418`) | ✗ |
| `MovingFrontierWindow` | `_CellSolve/_CellResidual → cell_kernel_batch` (`:887`/`:949`) | ✓ |
| `FullFieldWavefront` | same (`:1090`/`:1148`) | ✓ |

**What is and isn't single-sourced:** `cell_balance.py`'s denom/numer (`denom = 2|μ|·A_down +
(ΔA/w)·c_out + Σ_t·V`, additive **spatial + angular + collision**) IS single-sourced. NOT
single-sourced: (a) the **closure** `ψ_out=2ψ̄−ψ_in` / `ψ̄=½(in+out)` — written in ~10 sites
(`diamond.py:156/327/363`, `scan.py:331/332`, `loss_representation.py:1418/2004/2167`,
`operator.py:469/960`); (b) the **affine-scan form** `a=2|μ|A_total/denom−1` — written in 3
(`sweep_cache.py:419-456`, `ScanMarch._sweep_interior:1301`, `scan.py` crosswalk).

**It is the ERR-026 debt class.** Manifestation #7 was sweep-vs-matvec computing one operator via
two transcriptions that silently diverge on non-flat ψ. The ~10 closure sites agree only by
test-pinned coincidence (`rtol=1e-13` dual-view, `window≡full` oracle, nULP ScanMarch crosscheck);
the moment one closure site is "improved" (positivity fixup, WDD-τ, EC branch) without the others,
manifestation #7 reproduces **by construction**. The wavefront family is already L21-unified (both
route `cell_kernel_batch` → cannot drift); the scan family is the **unfinished half of the #196
unification**. The refactor makes the bug class **unspellable** (one closure site).

**Space↔angle is cleanly separable** (Phase 2.11): the M-M angular term enters `cell_balance` as a
closure-produced additive argument (`angular_denom_term`/`angular_numer_upstream`); slab's
`IdentityAngularClosure` returns zeros. The only fusion is the affine transform's `−1` (the closure
absorbed into the multiplier) — which is exactly what the new capability should own.

## The conceptual frame — schedule ⊗ cell-update (user, Nexus-confirmed)

The four strategies are the **sweep-schedule axis** (HOW cells are traversed/solved):
`CumprodScan` (Blelloch closed-form chain scan) / `ScanMarch` (x-scan ∘ y-march) / the
`_DAGWavefront` family (anti-diagonal independent-set walk). It is ALREADY a first-class selectable
mechanism (`LossRepresentation` Protocol + `default_for(mesh)` selection). The **cell-update axis**
(WHAT each cell computes: DD/LD/EC) is a separate selectable scheme (`CellUpdate` on `SNMesh`).
These should be ORTHOGONAL — any schedule × any cell-update — a (schedule ⊗ cell-update) tensor
product nested inside #236's space⊗angle theme.

**Confirmed (Nexus + instantiation-site grep):** the cell-update currently serves ONLY the
wavefront schedule. `cell_kernel_batch` is the *independent-batch solve* ("given cells whose
incoming faces are all known, solve them") — exactly the anti-diagonal's access pattern; only
`MovingFrontierWindow`/`FullFieldWavefront` (via `_CellSolve`/`_CellResidual` at `:887`/`:949`/
`:1090`/`:1148`) consume it. The scan/march schedules CANNOT use it: their cells are coupled along
the sweep direction (cell *i* inflow = cell *i−1* outflow), so they need the closed-form affine
recurrence `a=2sA/D−1` — which the cell-update does NOT expose, so they inline DD. (Nexus `callers`
on `cell_kernel_batch`/`_CellSolve.cell` returns empty — a dynamic-dispatch false negative via the
abstract `level_op.cell(...)` call; the instantiation-site grep is the definitive evidence.)

## The design — make the cell-update schedule-agnostic (TWO capabilities, one per schedule family)

The scan family (Blelloch parallel-reduce, `cumprod·(ψ₀+cumsum(b/cumprod))`) and the wavefront
family (per-level independent-batch solve) are **genuinely different folds** — they are the two
schedule families, and each needs a DIFFERENT thing from the cell-update. Collapsing them — forcing
the scan to call `cell_kernel_batch` per-cell — destroys the O(N/log N) advantage that is the
entire reason `CumprodScan`/`ScanMarch` exist. ⟹ the cell-update exposes ONE capability per schedule
family (both owning the SAME closure), so EVERY schedule consumes the scheme:

- KEEP `cell_kernel_batch` / `residual_kernel_batch` (wavefront family — already consumed).
- ADD to the `CellUpdate` Protocol (+ `CellUpdateBase` `NotImplementedError` defaults, mirroring the
  existing `cell_kernel_batch` default):
  - `is_affine_scannable: bool` (class trait — DD `True`, future Step/EC also affine; LD `False`).
  - `affine_scan_coefficients(...) -> AffineScanCoeffs` — returns the σ_t-epoch, **source-INDEPENDENT**
    pair `(a_attenuation, inverse_denom)`, where `a = 2|μ|·A_total/denom − 1`,
    `denom = 2|μ|·A_down + dA_w·c_out + Σ_t·V`. **Source of truth = the exact math currently in
    `CollisionCache.from_geometry` (`sweep_cache.py:419-453`)** — lift it verbatim into DD.
    **Inputs = the geometry arrays it reads** (`abs_mu, A_down, A_total, dA_w, c_out, V`) **+ `sig_t`**;
    pass the PRIMITIVE arrays, NOT the `GeometryCoefficients` type, so DD/`cell_update` stay decoupled
    from the `sweep_cache` layer (import direction). **Vectorized whole-tensor (the 3 numpy ops),
    NEVER per-cell** (per-cell would re-introduce the loop the L15 cache eliminated).
    ⚠ NOT the `cell_kernel_batch`-style `(s_axes, sigt_cells, Q_cells)` signature — that is the
    batch/wavefront capability; this is the scan capability. The per-sweep source term
    `b = 2·(QV + angular_contrib)·inverse_denom` is computed in the SWEEP BODY (`_run_1d_sweep`) from
    `inverse_denom` + Q — source-dependent, so NOT in this capability. `cumprod_a` stays in
    `CollisionCache` (the scan's transform of `a`, not DD math).
  - `scan_closure(face_in, face_out) -> psi_avg` — DD: `0.5*(in+out)` (the solve-direction inverse of
    the `2ψ̄−ψ_in` closure).
  - `transverse_close(psi_bar, psi_y_in) -> psi_y_out` — DD: `2*psi_bar − psi_y_in` (ScanMarch's
    y-march; a stub in PR1, wired in PR3).
- `DiamondDifference` implements all three (lifting the affine derivation from
  `sweep_cache.from_geometry` + the closures from the inline sites).
- A non-affine scheme (LD) sets `is_affine_scannable=False`, raises on `affine_scan_coefficients`;
  `ScanMarch.supports(mesh)` / `CumprodScan.supports` gain the `mesh.cell_update.is_affine_scannable`
  check; `default_for` falls through to `FullFieldWavefront` (the d-generic spine that needs only
  `cell_kernel_batch`). **Selection machinery unchanged** (existing `supports → Compatibility →
  default_for` fall-through, `:1450-1488`).

**Count-the-concepts test passes:** ~10 closure sites → 1 (DD); 3 affine transcriptions → 1 (DD);
+1 justified capability axis (scan vs batch). `cell_balance.py` keeps the additive denom; the L15
two-stratum cache keeps storage/lifetime (`CollisionCache.from_geometry` calls
`DD.affine_scan_coefficients`) — DD owns the math. This converts the coincidental
`test_cache_populator_matches_cell_balance_terms` (`rtol=1e-14` tripwire) into a true contract.

## Sequence — full uniform, slab-first (de-risk order)

This is NOT a "wait for ≥2 instances" case: the concept (DD's closure) is ESTABLISHED with 4+
consumers TODAY and the duplication is ACTIVELY re-spelling the ERR-026 class — the
`defer-abstraction-is-anti-no-benefit` exception applies (benefit = illegal-states-unrepresentable
for the twin-path bug class). Do the FULL uniform refactor, in de-risk order:

- **PR1 — slab (1-D Cartesian) + the matvec twin.** Add the DD capabilities; route
  `_run_1d_sweep:2004` AND `operator.py _compute_LpC:469` through them. Simplest path (no M-M);
  proves the capability against the bit-identity-critical slab regression snapshots (`rtol=1e-12`).
- **PR2 — curvilinear 1-D (`:2167`).** M-M is a separable additive arg (the closure keeps owning
  M-M). Highest correctness risk (ERR-058-adjacent) → second, after the capability is proven.
- **PR3 — ScanMarch 2-D (`:1301`/`:1418`) + `_compute_decomposition:960`.** Lowest perf risk
  (per-row recompute, no cache interaction; `#206`).

Each PR: DD bit-identical (or principled-equiv nULP where association order legitimately differs) +
the seam-honored gate goes green + qa + elegance-enforcer + ff-only merge. Main-agent writes inline
(surgical SN carve — NOT method-implementer, per standing guidance).

## Perf-neutrality + THE ONE TRAP

- Perf is preserved IFF `affine_scan_coefficients` stays the **vectorized cache builder** (same 3
  numpy ops, owned by DD instead of `sweep_cache`), never per-cell.
- ⚠ **THE TRAP (operation order):** `affine_scan_coefficients` must reproduce the CACHE's factoring
  (`b = 2·QV·inv_denom`, V folded into the collision term), NOT `cell_kernel_batch`'s explicit
  left-fold order — else slab snapshots break at `rtol=1e-12`. The single-sweep `--capture-baseline`
  bit-identity gate catches this; it is the chief correctness risk of PR1.

## Verification (test-architect spec — holds + extends to all strategies)

- **Single-sweep / single-matvec bit-identity** (`--capture-baseline`, `kind="direct"`, random ψ,
  `array_equal` rtol=0). End-to-end DD snapshots ALREADY drift on clean `main` (FP jitter, up to
  272k ULP) → bit-identity MUST be pinned at the single-sweep level, geometry-agnostic.
- **The probe → permanent "every strategy honors the scheme" gate** (now all 4, not just
  curvilinear). RED today → land WITH the fix (not a standing xfail). Keep +/− controls. vv Mode 8:
  `np.testing`/`pytest.fail` (fires under `-O`).
- **DD-via-new-capability ≡ DD-via-old-path** principled-equiv (nULP), ≥2g + heterogeneous +
  non-flat + non-zero inflow (anti-patterns #3/#4, Mode 9).
- **`test_cache_populator_matches_cell_balance_terms`** becomes a true round-trip contract.
- Ride existing pins: `test_diamond.py`, `test_streaming_operator_decomposition.py`,
  `test_wavefront_cumprod_equivalence.py`, curvilinear MMS + `test_keff_curvilinear.py`. Snapshot
  verdict: ZERO regen expected (any iterated snapshot past `SAFETY×conv_tol` = a bug).

## Downstream (closes the user's two gaps)

After the seam is uniform: **#158 LD/SC occupant** (`is_affine_scannable=False` → routes via
`FullFieldWavefront`) lifts the #233 pole-cell O(h)→O(h²); then demonstrate genuine independent
selection (DD-space ⊗ M-M-angle vs LD-space ⊗ M-M-angle). That fully closes gap (1) "the LD option"
and gap (2) "independent space⊗angle instantiation."

## Files
- `orpheus/sn/spatial/cell_update.py` (Protocol + Base capability defaults)
- `orpheus/sn/spatial/diamond.py` (DD impl — the single closure/affine owner)
- `orpheus/sn/spatial/sweep_cache.py` (`CollisionCache.from_geometry` → `DD.affine_scan_coefficients`)
- `orpheus/sn/loss_representation.py` (`:2004`,`:2167`,`:1301-1310`,`:1418`; `*.supports` gates)
- `orpheus/sn/spatial/scan.py` (closure sourced from scheme; `_scanmarch_row:331-332`)
- `orpheus/sn/operator.py` (matvec twin `:469`,`:960`)
- New tests (test-architect): `tests/sn/sweep/core/test_seam_honored.py`,
  `tests/sn/regression/test_sweep_seam_carve.py`
- Branches: `refactor/sn-cellupdate-seam-slab` (PR1) → `...-curvilinear` (PR2) → `...-scanmarch` (PR3)

## Design refinements (session-2 steering, code-confirmed)

- **`is_affine_scannable` is a CellUpdate trait, NOT a strategy bool.** There is no existing
  DAG-vs-DAG-free bool on the `LossRepresentation`s — the distinction is class-structural
  (`_DAGWavefront` subclass) and compatibility lives in each strategy's `supports(cls, mesh) ->
  Compatibility` (only mesh-geometry traits `is_1d`/`is_cartesian`/`ndim` today).
- **affine ≠ DAG-free — it's a compatibility relation.** A DAG-free (scan/march) schedule REQUIRES
  an affine cell-update; a DAG (wavefront) schedule works with ANY cell-update. An affine scheme
  (DD) runs on BOTH (CumprodScan production + FullFieldWavefront oracle = the sanctioned pairing).
  ⟹ `is_affine_scannable` is the PRECONDITION for selecting a scan schedule.
- **The check lives in the scan strategies' existing `supports`** (`CumprodScan.supports`/
  `ScanMarch.supports` add `and mesh.cell_update.is_affine_scannable`), enforced by the EXISTING
  `__post_init__` guard (`:312`) + `default_for` fall-through (`:1481`, first-supporting-strategy).
  No new machinery; illegal pairing (non-affine + scan) becomes unrepresentable.
- **CellUpdate is instantiated OUTSIDE the LossRepresentation** — on `SNMesh` (`cell_update=`,
  `geometry.py:271`); the strategy (frozen dataclass holding `mesh`) reads `self.mesh.cell_update`.
  So NO instantiation-time affine check; the check is at strategy SELECTION (the guard/`default_for`).
- ⚠ **Curvilinear + non-affine has NO strategy yet** (PR2 design point, not PR1): `default_for`
  falls a non-affine scheme through to the DAG family, which is Cartesian-only (`FullFieldWavefront`
  requires `is_cartesian`). So slab+non-affine → `FullFieldWavefront` ✓, but curvilinear+non-affine
  (the #233 LD payoff) → `IncompatibleRepresentation`. PR2 must add a **curvilinear non-affine
  schedule** — seed = the per-cell `cell_update.update` path at `_run_1d_sweep:2127` (the
  degenerate-ordinate vehicle). PR1 (slab) is unaffected.
- **Angular curvature term (`dA_w·c_out`) deferred to §3** — open: does any OTHER selectable
  curvilinear angular closure need this form, or is it M-M-specific? (PoleAngularClosure is its own
  selectable axis; M-M is one occupant.)

## DO NOT (study caveat)
Do NOT collapse the two access patterns (scan vs batch) into one capability chasing "elegance" —
they are different folds; forcing scan→`cell_kernel_batch` per-cell kills the Blelloch advantage.
Two capabilities on one scheme IS the elegant endpoint (the concept count shrinks).
