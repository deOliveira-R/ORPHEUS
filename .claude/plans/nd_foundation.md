# N-D General Foundation — dimension-generic SN spatial machinery (generalizes to 3D)

**Status:** FORWARD plan for a FUTURE session — picked up **AFTER Wave O
completes** (WavefrontFlux foundation A→B → SI Gauss-Seidel recovery → O.2).
Locked in 2026-06-04 at user request ahead of a context compaction. Details
firm up at pickup; this captures the goal, the structure, the generalization
targets, and the justification so a cleaned-context session resumes without
insight loss.
**Branch:** continue on `refactor/field-role-typing` lineage (or its successor).
**Mode:** main-agent direct authorship, turn-by-turn steering. Standing commit
exclusions as always (`.claude/agent-memory/`, `.mcp.json`, hooks,
`derivations/diagnostics/`).

---

## 0. The goal (one line)

**Tighten the SN spatial machinery into ONE dimension-generic foundation
(`d` axes) and validate it by adding 3-D Cartesian** — collapsing the
1-D + 2-D code paths into a single `d`-parametric implementation, so 3-D is a
*parameter* (`d=3`), not a new code path.

This is the `feedback_architecture_forward_not_legacy_fit` "extend to N-D, not
N+1" applied at the spatial-sweep layer, and it is *licensed now* by
`feedback_unify_after_two_instances`: **1-D and 2-D are the two working
instances; the d-generic unification is the third-instance step, with 3-D as
the validation it generalized** (not a premature guess).

---

## 1. Why it's tractable — the cochain frame is dimension-agnostic

(Established in `wavefront_flux_foundation.md` + the cross-domain-attacker memo
`.claude/agent-memory/cross-domain-attacker/field_role_typing_faceflux_frames.md`.)
Every spatial concept generalizes by "add an axis", because none has `d` baked
into the *math*:

| concept | 1-D | 2-D | 3-D | generic form |
|---|---|---|---|---|
| cell-average ψ̄ (d-cochain) | per cell | per cell | per cell | dimension-free |
| face flux `WavefrontFlux` (codim-1 cochain) | 2 pts/cell | 4 edges/cell | 6 faces/cell | per-axis: `d` orientations |
| DD closure `ψ_out=2ψ̄−ψ_in` | x | x,y | x,y,z | `Σ_axis sx_axis·ψ_in_axis` |
| upwind DAG level | `i=k` | `i+j=k` | `i+j+k=ℓ` | `Σ_axis idx_axis = ℓ` (diagonal hyperplane) |
| octants (2^d) | 2 | 4 | 8 | sign-tuple over `d` axes |
| boundary trace `ι*` (2d faces) | 2 | 4 | 6 | codim-1 edge subset, 2 per axis |
| streaming coeff `2|μ_axis|/Δaxis` | 1 | 2 | 3 | per-axis vector |

**Precondition (delivered by the WavefrontFlux foundation §3a′):** the
`WavefrontFlux` / interior `FaceLayout` / `ι*` **types are already built
axis-parametric** (face field over `axes`, layout keyed by `(axis, position)`,
`ι*` over each axis's edge subset), with a unit test proving the API runs
`axes=(0,)` and `axes=(0,1)` through one code path. So this session generalizes
the *call sites* and *strategies*, not the field types.

---

## 2. Generalization targets (the 2-D/1-D-hardcoded sites to make `d`-generic)

All `file:line` are at the WavefrontFlux-foundation HEAD; re-confirm at pickup
(they will have moved). These are the sites that still hardcode "2 axes / 4
faces / 4 octants / explicit x,y":

1. **`SweepDependencyGraph.from_cartesian_2d`** (`sweep_graph.py:191`, levels
   `i+j=k`) → **`from_cartesian(d)`**: the upwind DAG over the `d`-dimensional
   index lattice (levels `Σ idx = ℓ`), `2^d` octant graphs, per-octant index
   reversal generalized over the sign-tuple. `_make_slice` generalized to
   `d`-axis face-index wiring.
2. **`DiamondDifference.update_batch` / `residual_batch`** (`diamond.py:273/361`,
   explicit `psi_x`/`psi_y`, `sx·ψ_in_x + sy·ψ_in_y`) → **sum over `d` axes**:
   `denom = Σt + Σ_axis sx_axis`, `ψ̄ = (Q + Σ_axis sx_axis·ψ_in_axis)/denom`,
   the `d` face-pairs per cell. The closure `ψ_out = 2ψ̄ − ψ_in` is already
   per-axis — vectorize over `axes`.
3. **`_sweep_2d_wavefront`** + **`_run_1d_sweep`** → **ONE `_sweep_wavefront`**
   over `d` (the tightening: the 1-D and 2-D sweep bodies COLLAPSE into the
   generic walk — the unify-after-two payoff; retire the two specialized
   bodies). The matvec twin `_apply_2d_cartesian` + the 1-D
   `_compute_LpC`/`_compute_decomposition` likewise collapse onto the
   `d`-generic `graph.residual`.

   **⚠ The 1-D-fold seam (evaluated 2026-06-04 — explorer memo
   `wavefront_flux_1d_tightening_verdict.md`).** This collapse is the ONLY
   place the WavefrontFlux infrastructure reaches 1-D — and it is NOT a free
   retrofit, because **1-D today is a parallel-prefix SCAN, not a wavefront**:
   - 1-D interior face fluxes = `ordinate_scan`'s transient `(nx, K, ng)`
     **chain-ordered, K-batched** output (`sweep.py` slab path) — the
     structural antithesis of the cochain's `(N, ng, nx+1)` (cells trailing,
     un-chained, ordinate-leading). The 1-D matvec (`_compute_LpC`) has no
     interior buffer at all (a per-cell stack local threaded through the DAG
     walk). So `WavefrontFlux` does NOT fit the 1-D scan as-is (confirmed
     WRONG-FIT — typing it would fight the L16 cumprod for zero gain).
   - **What unification harvests** (only AFTER both folds become one
     `dag_walk`-driven traversal): (a) `WavefrontFlux` + `seed`/`absorb`
     (ι_*/ι*) become 1-D's natural interior representation too — a
     re-expression, not a literal swap; (b) the **boundary-trace exchange**
     (the one genuinely-shared concept — both folds seed inflow / persist
     outflow on the same `BoundaryFlux`, but 1-D does it per-direction-chain
     via a scan-tail `[...][-1][ords]` while 2-D does `wavefront.seed/absorb`)
     unifies into ONE typed `TraceExchange`; (c) the DD closure
     (`0.5·(ψ_in+ψ_out)` scan vs `2ψ̄−ψ_in` per-level) collapses to one
     `C¹→C⁰` averaging primitive.
   - **HARD CONSTRAINT on the collapse:** the d=1 parallel-prefix scan is the
     *better* algorithm for the linear recurrence (Blelloch closed form,
     "2 scan calls per sweep"), NOT legacy. A naive `d`-generic
     forward-substitution walk is O(nx) **sequential** for d=1 — a regression.
     So the unified `_sweep_wavefront` MUST either keep the scan realization
     for d=1 (the walk specializes its fold per `d`) OR express the
     `d`-generic walk itself as a segmented/parallel scan that recovers the
     d=1 cumprod. Pin: 1-D wall-clock must not regress (L16) through the
     collapse. This is the load-bearing design question of the collapse.
4. **`Mesh2D`** (`geometry/mesh.py:491`) + `Mesh1D` → a `d`-generic Cartesian
   mesh (or `Mesh3D` as the third instance, then unify). `SNMesh` streaming
   coefficients over `d` axes (`str_x`, `str_y` → `str[axis]`).
5. **Quadrature** (`directional.py`): 3-D needs the FULL sphere (8 octants,
   level-symmetric / product over `(μ,η,ξ)` with `μ²+η²+ξ²=1`); 2-D used the
   half-sphere reduction. Generalize `quad.octants` + `reflection_index(axis)`
   over `d` axes. (This is the heaviest non-spatial piece — the angular
   discretization in 3-D.)
6. **Boundary**: `BoundaryFlux` / `FaceLayout` over `2d` faces
   (`+zmin/zmax`) — the same `FaceLayout` generalization the interior needs;
   `SNBoundaryOperator` over `2d` faces (already loops `layout.faces`, so mostly
   free once the layout is `d`-generic).

   **⚠ Mesh-side caveat (issue #220).** "Mostly free" covers only the
   `FaceLayout` storage + the `SNBoundaryOperator` loop. The MESH-side BC
   *resolution* surface is NOT yet `d`-generic: `SNMesh._resolve_bcs`
   (`orpheus/sn/geometry.py`) is a `Mesh1D | Mesh2D` `isinstance` split that sets
   named `bc_xmin`/`bc_xmax`/`bc_ymin`/`bc_ymax` (+ `bc_left`/`bc_right`)
   attributes — which 3-D would force to `bc_zmin`/`bc_zmax`. Replace with a
   `FaceLabel`-keyed `SNMesh.bc` dict + a single `face_labels` loop (R-1 C4's
   unfinished half) BEFORE or AS part of the 3-D boundary step. See issue #220.

---

## 3. Scope + non-scope
- **In scope:** 3-D **Cartesian** SN (the natural extension of the 2-D
  wavefront). The `d`-generic unification of 1-D + 2-D + 3-D Cartesian.
- **Out of scope (separate future beasts):** 3-D **curvilinear** (spherical /
  cylindrical in 3-D — the curvilinear closures are currently 1-D only; the
  Carlson/M-M pole machinery is a different generalization axis). MoC stays on
  its fiber-bundle/sheaf path (`project_moc_structure` — NOT this wavefront DAG).
- **Storage:** by this session, WavefrontFlux path **(B) (moving-frontier
  window)** has landed — and 3-D is exactly where it pays most (the
  volume↔surface ratio makes the full-field `O(N·ng·nx·ny·nz)` painful; the
  frontier window `O(N·ng·(nx·ny+ny·nz+nx·nz))` is the win). So B is a
  *prerequisite* for practical 3-D, not just an optimization.

---

## 4. Verification shape (firm up at pickup)
- **3-D references:** a 3-D Cartesian MMS (manufactured solution; extend
  `derivations/continuous/mms/sn.py`); 3-D homogeneous reflective k_inf (= the
  same `νΣ_f/Σ_a` anchor, dimension-free → must still be 1.875 for mixture A 2g);
  3-D per-ordinate streaming-equilibrium flat-flux (the L0 anchor, dimension-free).
- **Regression (the tightening must not break the existing instances):** the
  `d`-generic code MUST reproduce 1-D + 2-D **bit-identically** (the unification
  is a refactor, not a numerics change — same gates: octant snapshots, Gate-K,
  slab/curvilinear MMS, sentinel). This is the load-bearing constraint — collapse
  the code paths WITHOUT moving a bit on 1-D/2-D.
- **Perf (L16):** the `d`-generic walk MUST stay vectorized (the `Σ_axis` is a
  loop over `d≤3`, not over cells/faces); full-suite wall-clock not regressed.
- Default `python -O`; NO `continuous_get` (#212). Bit-identity discipline per
  `vv-principles`.

## 5. Risk ranking
1. **The tightening regresses 1-D/2-D (HIGH).** Collapsing two proven sweep
   bodies into one `d`-generic walk risks a subtle 1-D/2-D bit-move. Mitigation:
   bit-identity gates on BOTH existing instances at every step; the `d`-generic
   walk is validated against the retiring specialized bodies before they're
   deleted (`feedback_retirement_means_test_migration`).
2. **Vectorization loss (HIGH, L16).** A `d`-generic walk that loops over
   cells/faces instead of `axes` → catastrophic regression. The `Σ_axis` must be
   a tiny `d≤3` loop over vectorized per-axis arrays.
3. **3-D quadrature (MEDIUM-HIGH).** The full-sphere 8-octant angular set is the
   heaviest new piece; may need a `literature-researcher` + `test-architect` pass
   on the 3-D level-symmetric / product quadrature + its `reflection_index`.
4. **Memory without (B) (MEDIUM).** If (B) the frontier window isn't in by then,
   3-D full-field memory may be prohibitive — confirm (B) landed first.
5. **Curvilinear entanglement (LOW if scoped out).** Keep 3-D Cartesian clean of
   the 1-D curvilinear pole machinery; don't let the `d`-generic refactor
   accidentally couple them.

## 6. Pickup checklist (fresh session, AFTER Wave O)
1. Confirm Wave O complete: WavefrontFlux (A) + (B) landed, G-S recovery landed,
   O.2 landed (the loss-operator driver + adjoint). The `WavefrontFlux` /
   `FaceLayout` / `ι*` types are axis-parametric (the §3a′ unit test green).
2. Read this plan §1–§2 (the generic forms + the generalization targets) and the
   cross-domain-attacker cochain memo. Dispatch **explorer** to re-confirm the
   §2 site line numbers (they will have moved) + **test-architect** for the
   3-D reference + 1-D/2-D bit-identity regression plan (proactive: a carve that
   collapses code paths across dimensions).
3. Tighten FIRST (unify 1-D + 2-D onto the `d`-generic walk, bit-identical, retire
   the specialized bodies), THEN extend (add 3-D Cartesian: `Mesh3D`, 8-octant
   quadrature, the 6-face boundary, the 3-D references). Turn-by-turn; per-step
   bit-identity + wall-clock gates; elegance-enforcer per phase.
4. Docs (archivist): the `d`-generic cochain foundation in the theory pages;
   retire the dimension-specific narrative; the storage×role grid (#205) FACE
   locus now `d`-generic. Update `wave_o_operator_typing.md` lineage + auto-memory.

## 7. References
- `.claude/plans/wavefront_flux_foundation.md` (§3a′ — the axis-parametric types
  this session builds on; the cochain frame).
- `.claude/agent-memory/cross-domain-attacker/field_role_typing_faceflux_frames.md`
  (the dimension-agnostic cochain detection + the `(octant×face)` graph).
- `.claude/plans/si_gauss_seidel_recovery.md` (the G-S schedule, which must
  also be `d`-generic — the `(octant×face)` graph is already dimension-free).
- `feedback_architecture_forward_not_legacy_fit`, `feedback_unify_after_two_instances`
  (the governing principles — and the tension they resolve: two instances now
  license the N-D unification, 3-D validates it).
- `feedback_retirement_means_test_migration` (retire the 1-D/2-D specialized
  sweep bodies AS the `d`-generic walk replaces them).
