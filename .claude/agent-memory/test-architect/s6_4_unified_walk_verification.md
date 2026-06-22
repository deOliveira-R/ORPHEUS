---
name: s6-4-unified-walk-verification
description: S6.4 V&V gate plan — the unified kernel-parameterized octant walk (_OctantWalk, d-generic) + DAG-ownership move (#222). Sweep + matvec walk ONE octant frame, forked only at a CELL-KERNEL object (apply_windowed/residual_windowed) + an EMIT policy, NOT a boolean. Anchor set (window≡full sweep+matvec, scan-march G2.c, affine-carve golden, A2D-1, SI≡Krylov≡k_inf, loss_action convention pin). The one-walk discriminating spy test. A2D-1 RETIRE-in-favour-of-output-oracle recommendation. DAG-ownership: mesh-no-longer-exposes-sweep_graphs structural test + geometry.py:343/350 None-slot gone. Mode-8/Mode-9. Per-sub-step staging (a)matvec-share→(b)sweep-into-walk→(c)DAG-move→(d)fold-full-field. ADDENDUM 2026-06-11: sub-steps (e) collapse sweep_graph's 4 level-walks→2 direction-parameterized walks (CellUpdate override point becomes the pure kernel pair; per-level-emit-reorder + left-fold-order bit-id risk) + (f) module geography — sweep.py DISSOLVES (scan→spatial/scan.py, walk+driver+transport_sweep→loss_representation.py, SweepDependencyGraph.for_shape classmethod, lazy-import cycle gone, ~30 test files migrate, Sphinx :func: role migration). Sequencing: (c) should ALREADY mint for_shape on the graph class so (f) doesn't move it twice. Extends [[s6-relayering-verification]].
metadata:
  type: project
---

# S6.4 — Unified kernel-parameterized octant walk + DAG ownership (V&V gate plan)

S6.4 is the DEEPEST L21 realization of the SweepStrategy carve (#222): `sweep`
and `loss_action` (matvec) walk the SAME octant DAG / SAME frontier / SAME
boundary seed-shed — the face recurrence `out = 2ψ̄ − in` is LITERALLY identical
— forking ONLY at (i) the **cell kernel** (`cell_kernel_batch` SOLVE vs
`residual_kernel_batch` APPLY) and (ii) the **emit policy** (sweep → ψ +
scalar/moment projection; matvec → residual bulk + O.4b boundary defect). Today
the octant orchestration is written **multiple times in lockstep**; S6.4 extracts
ONE `_OctantWalk2D` shared by sweep + matvec, and moves the DAG (`sweep_graphs`)
off the mesh into `_DAGWavefront`. Proactive `test-architect` gate (lesson L17)
BEFORE any code moves — the carve crosses sweep↔matvec↔mesh↔representation.

Extends [[s6-relayering-verification]] (the S6 anchor vocabulary, the A2D-1
regenerate-in-commit protocol, Mode-8/Mode-9 discipline) and
[[sweep-strategy-carve-verification]] / [[scan-march-verification]]. Worktree
`worktree-sn-nd-layout`, HEAD `7bbeb5f` (S6.3b LANDED — the walk is already off
the operator into `loss_representation.py`; S6.4 collapses the Fork-B1 lockstep
twins INSIDE the representation module). Authoritative design =
`.claude/plans/sn_sweep_strategy.md` §S6.4 (the two "S6.4" bullets, committed
`7bbeb5f`).

## EXECUTIVE SUMMARY (5 lines)
1. **Anchor set (all VERIFIED EXISTS @ 7bbeb5f, collect `-O`-clean, 39 tests):**
   `window≡full` SWEEP **and** MATVEC oracle (`test_2d_full_field_oracle.py`, the
   load-bearing both-directions proof the shared walk reproduces sweep AND matvec),
   sweep-graph `window≡full` d1/d2/d3 (`test_sweep_graph_window_equivalence.py`),
   ScanMarch G2.c (`test_scan_march_equivalence.py`), affine-carve golden
   (`test_affine_carve_bit_identity.py`), `SI≡Krylov≡k_inf` (`test_keff_2d.py`),
   `loss_action` convention pin (`test_loss_action_convention.py`), A2D-1 source-hash.
2. **The one-walk discriminating test (the structural payoff):** a SPY on the
   shared `_OctantWalk2D._interior_walk` hit by BOTH a `sweep` and a `loss_action`
   call → assert SAME code object exercised in both directions. FAILS today (3×
   duplicated octant frames → no shared walk), PASSES after S6.4. `-O`-safe
   (`pytest.fail`). Stage `xfail(strict=True)`→xpass.
3. **A2D-1: RETIRE in favour of the output oracle (recommendation).** A2D-1 pins
   `sha256(getsource(MovingFrontierWindow.loss_action))`. S6.4 moves that body's
   octant frame + boundary block INTO `_OctantWalk2D` → the hash regenerates. Once
   the walk is SHARED and covered by the `window≡full` MATVEC output oracle, the
   source-hash tripwire is REDUNDANT (its job — "did someone touch this body" —
   transfers to the bit-identical output oracle that survives the relocation). RETIRE
   it at S6.4 (delete the class), document the handoff. (Fallback if kept: retarget
   to `_OctantWalk2D._interior_walk` + regenerate in-commit.)
4. **DAG-ownership structural test:** assert the mesh NO LONGER exposes
   `sweep_graphs` (`not hasattr(SNMesh-instance, "sweep_graphs")` for ALL coord
   systems) AND `CumprodScan`/`ScanMarch` never reference it (grep-gate +
   construction-gate) — illegal-states-unrepresentable; plus the `geometry.py:343/350`
   curvilinear `self.sweep_graphs = None` None-slot is GONE; plus per-mesh-shape DAG
   cache produces byte-identical graphs (the window≡full sweep+matvec stays green
   after the move).
5. **Staging (each independently bit-identity-gated):** (a) extract `_OctantWalk2D`
   for the two MATVEC variants (window + scanmarch) first — cheapest, A2D-1 fires;
   (b) bring `sweep_octant_group` + `_sweep_2d_scanmarch` into it (the deepest L21
   unification, the one-walk test flips); (c) the DAG-ownership move (mesh → family);
   (d) fold `FullFieldWavefront` (the third octant frame + boundary block). Per-stage
   gate table = §7.

## Claim-layer / pillar gate (vv-principles §1.5 — MUST pass before writing)
Every S6.4 gate is a **behavior-PRESERVING carve gate**, NOT a new-physics claim —
the converged answer does not change (a re-org of orchestration, S6.4 forks no
default). Claim layers + pillars:
- **Bit-identity anchors** (window≡full sweep+matvec, sweep-graph window≡full,
  affine-carve golden, A2D-1 while kept): NOT a pillar claim — *free verification
  by inheritance* (`np.array_equal`/sha256, exact).
- **ScanMarch G2.c** (`ScanMarch.loss_action ≡ FullFieldWavefront.loss_action`):
  flux-shape claim → **`FullFieldWavefront` oracle** (a different schedule of the
  SAME lower-triangular solve), transitively pinned to k_inf (L8) + φ=Q/Σ_t (L7).
  Principled-equiv (nulp) by construction (row-march vs anti-diagonal differ at
  FP-association). Structurally independent. PASS.
- **`SI≡Krylov≡k_inf`** (G6): EIGENVALUE claim → **closed-form pillar**
  (k_inf=λ_max(A⁻¹F) transfer matrix). ≥2G non-flat (cardinal rule — NO 1G
  eigenvalue test). MMS NOT used for the eigenvalue leg. PASS.
- **Curvilinear angular-transpose** (rides §5 of [[s6-relayering-verification]]):
  reciprocity `⟨Aψ,φ⟩_G=⟨ψ,A.Hφ⟩_G` → **closed-form pillar**, independent `_g_inner`
  metric. PASS. (S6.4 touches the 2-D Cartesian octant frame ONLY — the curvilinear
  1-D transpose path `CumprodScan.loss_action_transpose` + `closure.angular_adjoint`
  is NOT re-carved here; `test_g_adjoint_reciprocity` stays as the regression sentinel,
  already Mode-8-fixed at S6.3a `ee91db4`.)
Chain of trust terminates in: closed-form k_inf transfer matrix, the
`FullFieldWavefront` d-generic spine (different schedule, same operator), the
analytical φ=Q/Σ_t + k_inf=1.875 grounds. NONE is another procedurally-different
derivation of the same identity. **Gate PASSES. Matrix ready to write.**

---

## 0. WHAT S6.4 IS (the carve the gates target) — and what stays invariant

S6.3b LANDED the walk OFF the operator INTO `loss_representation.py` (each
representation's `loss_action` OWNS its `(L+C)ψ` walk). But the octant
orchestration is now duplicated MULTIPLE TIMES IN LOCKSTEP — the documented Fork-B1
IOU. VERIFIED this session, the duplication inventory:

**MATVEC side (3 octant frames, in `loss_representation.py`):**
- `MovingFrontierWindow.loss_action` (`:451`) — octant projection + pure-z + face
  strings + `graph.residual_windowed` (`:534`) + the O.4b boundary block (`:547-558`).
- `ScanMarch.loss_action` (`:779`) — IDENTICAL octant projection + pure-z + face
  strings + the row-march `_x_scan_faces` interior walk + the O.4b boundary block
  (`:897-909`, the docstring SAYS "IDENTICAL to MovingFrontierWindow.loss_action").
- `FullFieldWavefront.loss_action` (`:611`) — d-generic octant projection + pure-z +
  `graph.residual` (full cochain) + the O.4b boundary block (`:689-699`).

**SWEEP side (matching octant frames, in `sweep.py`):**
- `sweep_octant_group` (`:796`) — octant projection + pure-z + face strings + inflow
  read + `graph.apply_windowed` (`:936`) + outflow shed into `boundary_flux`.
- `_sweep_2d_scanmarch` (`:1276`) — docstring (`:1323-1335`) SAYS its octant
  projection + boundary-trace I/O is "a DELIBERATE Pattern-2 duplication of
  `sweep_octant_group`'s scaffold — **edit both in lockstep**".
- `_sweep_full_field` (`:1451`) — the oracle sweep counterpart.

The face recurrence (`out = 2ψ̄ − in`) is the same; only (i) the **cell kernel**
(`cell_kernel_batch` SOLVE `diamond.py:273` divide-by-`D` with source `q`, ψ̄ UNKNOWN
| `residual_kernel_batch` APPLY `diamond.py:331` pure-reflection `α=−1` no divide,
ψ̄ GIVEN) and (ii) the **emit policy** (sweep → ψ angular/moment write into the field
| matvec → `LpC` bulk + the O.4b `streamed−given` boundary defect) differ.

**THE CARVE: extract `_OctantWalk2D` — the ONE octant traversal** (octant projection
+ pure-z `Q/Σ_t` branch + effective-sign + in/out face strings + frontier walk +
boundary seed/shed + the O.4b boundary block), **parameterized by a CELL-KERNEL
object** (`apply_windowed` | `residual_windowed` | the row-march scan) **+ an EMIT
policy** — **NOT a boolean `is_solve` flag** (the constraint, §6). `sweep` and
`loss_action` STAY the two NAMED public faces; both delegate to the ONE walk.

**The three "edit both in lockstep" IOU markers S6.4 retires** (retirement audit):
`loss_representation.py:479` (window↔scanmarch matvec), `loss_representation.py:808`
(scanmarch↔window matvec), `sweep.py:1328` (scanmarch-sweep↔sweep_octant_group).

**INVARIANT:** S6.4 changes NO converged answer and FORKS NO default
(`MovingFrontierWindow` stays the d=2 default, `CumprodScan` stays d=1, scan-march
stays opt-in). So bit-identity holds end-to-end (affine-carve golden byte-identical;
window≡full sweep+matvec `array_equal`; ScanMarch G2.c stays nulp-principled).

---

## 1. Bit-identity anchors that pin the unified walk as behavior-preserving

All paths VERIFIED to exist + collect `-O`-clean @ HEAD `7bbeb5f` (39 tests, this
session). Class: **bit-identity** (free inheritance) | **value-ground**
(structurally-independent reference) | **principled-equiv** (nulp).

| Anchor | File::test (worktree-relative) | Class | Gates sub-step | Why load-bearing for S6.4 |
| ------ | ------------------------------ | ----- | -------------- | ------------------------- |
| **window≡full SWEEP** | `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py::test_sweep_window_equals_full_field_end_to_end` (4 cases incl. NON-SQUARE 12×7, 5×9; `np.testing.assert_array_equal` of angular + scalar + moment) | bit-identity | **(b),(c),(d)** | **THE proof the kernel-parameterized walk reproduces the SWEEP direction.** After `sweep_octant_group`'s frame moves into `_OctantWalk2D`, this must stay byte-identical. NON-SQUARE catches an x↔y swap in the shared walk (a SQUARE octant snapshot is x↔y-blind, [[phase5a_5c_angular_windowing_lessons]] THE TRAP analogue). |
| **window≡full MATVEC** | `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py::test_matvec_window_equals_full_field_end_to_end` (same 4 cases; bulk residual + boundary-block residual, `assert_array_equal`) | bit-identity | **(a),(c),(d)** | **THE proof the shared walk reproduces the MATVEC direction.** After the window/full-field matvec octant frames move into `_OctantWalk2D`, byte-identical. The output-identity oracle the A2D-1 regen leans on (§2). |
| **sweep-graph window≡full d1/d2/d3** | `tests/sn/sweep/core/test_sweep_graph_window_equivalence.py::test_solve_window_equals_full_field` + `::test_residual_window_equals_full_field` (shapes incl. NON-SQUARE (16,64),(12,7),(5,9) + synthetic d=3 (3,2,3),(4,3,2)) | bit-identity | **(c),(d)** | d-generic window≡full at the GRAPH layer — pins the DAG-cache move (§3) produces byte-identical graphs. The d=3 synthetic admission proves `_OctantWalk2D` ⊕ DAG-cache stays general (no d=2 hardcode crept in). |
| **ScanMarch G2.c** | `tests/sn/sweep/cartesian_2d/test_scan_march_equivalence.py::test_scanmarch_residual_equals_oracle` + `::test_scanmarch_sweep_equals_oracle` + `::test_scanmarch_moment_equals_window` (`_NULP_BOUND`, NON-SQUARE 12×7/5×9, LS-4, 2G het aniso) | principled-equiv (nulp, cross-schedule) | **(a),(b),(c),(d)** | `ScanMarch.{loss_action,sweep} ≡ FullFieldWavefront` at nulp. When the scanmarch octant frame + boundary block move into `_OctantWalk2D` (the row-march interior walk stays the scanmarch `_interior_walk`), this stays nulp-green. NOT bit-id (FP-association by construction). |
| **affine-carve golden** | `tests/sn/solve/test_affine_carve_bit_identity.py::test_converged_flux_bit_identical_after_affine_carve` (3 params `si_2d_p1_aniso_het`/`krylov_2d_p1_aniso_het`/`si_slab_2g_het`; sha256 of converged `bulk.values`+`phi`; `-O`-safe `raise`) | bit-identity (end-to-end output) | **ALL** | End-to-end converged-flux byte-identity. S6.4 forks no default → byte-identical. If ever flips a default, REGENERATE per the Fork-B2 discipline ([[scan-march-verification]] §G5). |
| **SI≡Krylov≡k_inf** | `tests/sn/eigenvalue/test_keff_2d.py::test_si_krylov_heterogeneous_2g_nonflat_flux` + `::test_default_entry_hits_kinf` + `::test_2g_eigenvector` + `::test_homogeneous_exact` | value-ground (closed-form k_inf transfer-matrix; ≥2G non-flat) | **ALL** (the eigenvalue ground; both sweep + matvec converge to k_inf) | The structurally-independent VALUE the bit-identity oracles inherit from (vv §1.5: ULP-distance necessary-never-sufficient). |
| **loss_action convention pin** | `tests/sn/operators/test_loss_action_convention.py::test_loss_action_is_full_loss_LpC_flat_reflective` + `::test_apply_equals_loss_action_minus_independent_collision_het` (`[slab_2g]`,`[cart2d_2g]`) | value-ground + bit-identity glue | **(a),(b)** | Pins `loss_action` returns `(L+C)ψ` (NOT `L·ψ`) and `op.apply = loss_action − σ_t·ψ.bulk`. The shared walk MUST preserve this convention — `cart2d_2g` drives the 2-D octant frame the carve touches. |
| **L7 φ=Q/Σ_t d=2** | `tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py::test_2d_octant_sweep_closed_form_anchor` | value-ground | underpins G2.c transitively | STAYS. |
| **L8 k_inf=1.875 d=1** | `tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py::test_cumprod_path_hits_analytical_kinf` | value-ground | underpins G2.c transitively | STAYS. |
| **A2D-1 source-hash** | `tests/sn/operators/test_streaming_operator.py::TestT4dApply2DCartesianSourceHashPin::test_apply_2d_cartesian_source_hash_unchanged` (`sha256(getsource(MovingFrontierWindow.loss_action))`, `:1321/:1333`) | bit-identity (source text) | **(a)** (body moves) | **REGENERATES at (a)** — see §2 + the RETIRE recommendation. |

The two `window≡full` oracles (SWEEP + MATVEC) are the LOAD-BEARING pair: they prove
the ONE kernel-parameterized walk reproduces BOTH directions byte-for-byte. Every
S6.4 sub-step that touches the 2-D octant frame is gated by at least one of them.

---

## 2. The A2D-1 retarget — and the RETIRE recommendation

**Current state (verified):** A2D-1 pins
`sha256(inspect.getsource(MovingFrontierWindow.loss_action))`
(`test_streaming_operator.py:1321`+`:1333`; retargeted from `_apply_2d_cartesian` at
S6.3). S6.4 moves `MovingFrontierWindow.loss_action`'s octant frame + boundary block
INTO `_OctantWalk2D` — the body shrinks to a thin delegation
(`return self._octant_walk(kernel=residual_windowed, emit=...)`). The source text
changes → the hash MUST change.

**RECOMMENDATION: RETIRE A2D-1 at sub-step (a), in favour of the `window≡full`
MATVEC output oracle.** Rationale (the aggressive-retirement discriminator,
[[feedback_aggressive_retirement]]):

- A2D-1's JOB is "did someone touch this matvec body without thinking" — a tripwire
  on a body that was (pre-S6.3) procedurally-transcribed and NOT covered by a tight
  output oracle. After S6.4 the body is the SHARED `_OctantWalk2D._interior_walk`,
  and `test_matvec_window_equals_full_field_end_to_end` (`assert_array_equal`, 4
  cases incl. non-square) is a TIGHTER tripwire on the SAME concern: any drift in the
  shared walk's matvec output fails `array_equal` against the `FullFieldWavefront`
  oracle. The source-hash adds NOTHING the output oracle doesn't catch, and a
  source-hash on a SHARED walk is WORSE than on a private one — every legitimate
  refactor of the shared frame (which now serves sweep AND matvec AND scanmarch)
  trips a hash the author must hand-bump, with no behavior signal.
- This is NOT the "fuller-view oracle" exception (that keeps a genuine STRUCTURAL
  reference). A2D-1 is "same-math-pinned-by-text, structurally-covered-elsewhere" —
  exactly the retire case. The output oracle IS the structural reference.

**Retirement protocol:** in the sub-step (a) commit, DELETE
`TestT4dApply2DCartesianSourceHashPin`; add a one-line note in the
`test_2d_full_field_oracle.py` MATVEC test docstring recording it as the A2D-1
successor ("supersedes the retired A2D-1 source-hash — the output-identity tripwire
on the 2-D matvec body, now the shared `_OctantWalk2D` walk"). Cite the MATVEC oracle
green in the commit body as the output-identity evidence. Run the elegance-enforcer
retirement audit (§7) — A2D-1 gone, no dangling `getsource` ref.

**Fallback if the implementer/qa prefers to KEEP A2D-1** (less preferred): RETARGET
`inspect.getsource(MovingFrontierWindow.loss_action)` →
`inspect.getsource(_OctantWalk2D._interior_walk)` (the new shared home), regenerate
`EXPECTED_SHA256`, add a history line ("S6.4: matvec body relocated into the shared
`_OctantWalk2D._interior_walk`"), bump in the SAME commit as the body move. NEVER
pin the OLD `MovingFrontierWindow.loss_action` text after the body moves — the
retirement audit (§7 sub-step (a)) catches a stale pin: grep
`MovingFrontierWindow.loss_action` after the move; if the octant frame is gone, the
A2D-1 target MUST have moved or A2D-1 retired. (Per [[s6-relayering-verification]]
§2 the stale-pin failure mode is impossible-to-ship-silently — hash mismatch OR
AttributeError fires loudly; but the RETIRE path removes the brittleness entirely.)

---

## 3. The DAG-ownership verification

S6.4 moves the DAG (`sweep_graphs`) OFF the mesh INTO the `_DAGWavefront` family
{`MovingFrontierWindow`, `FullFieldWavefront`}, lazily built + cached per
mesh-SHAPE by an accessor the family OWNS. VERIFIED only `_DAGWavefront` uses it
(`sweep.py:910/1533` + `loss_representation.py:524/672`); the curvilinear
`geometry.py:343/350` `self.sweep_graphs = None` is the illegal-state None-slot.

### (a) Structural test — the mesh no longer exposes `sweep_graphs` (illegal-states-unrepresentable)

**File:** NEW `tests/sn/geometry/test_dag_ownership.py` (foundation; `-O`-safe).

```python
@pytest.mark.foundation
@pytest.mark.parametrize("coord", ["slab", "cart2d", "sphere", "cyl"])
def test_mesh_does_not_expose_sweep_graphs(coord):
    """[L0 structural] The DAG is OWNED by _DAGWavefront, NOT the mesh.

    Post-S6.4 the mesh is pure geometry — `sweep_graphs` lives on the
    representation family (cached per mesh-shape). A mesh attribute would be
    the illegal-state None-slot (geometry.py:343/350 curvilinear `= None`).
    """
    sn = _build_mesh(coord)                 # the standard per-coord builder
    if hasattr(sn, "sweep_graphs"):
        pytest.fail(
            f"{coord}: SNMesh still exposes `sweep_graphs` — the DAG must be "
            "owned by _DAGWavefront (cached per mesh-shape), not the mesh "
            "(geometry.py:343/350 None-slot is the illegal-state smell S6.4 closes)."
        )


@pytest.mark.foundation
def test_dag_free_reps_never_reference_sweep_graphs():
    """[L0 structural] CumprodScan / ScanMarch never touch the DAG.

    DAG-free representations are constructed on curvilinear / 1-D / row-march
    meshes where no anti-hyperplane DAG exists. They must NOT reference
    `sweep_graphs`/`OctantLabel` (would re-introduce the substrate coupling).
    """
    import inspect
    for rep_cls in (CumprodScan, ScanMarch):
        for name in ("sweep", "loss_action", "loss_action_transpose"):
            src = inspect.getsource(getattr(rep_cls, name))
            if "sweep_graphs" in src:
                pytest.fail(
                    f"{rep_cls.__name__}.{name} references `sweep_graphs` — a "
                    "DAG-free representation must never mention the DAG substrate."
                )
```

Note: the grep-gate (`inspect.getsource` for `sweep_graphs`) is `-O`-safe and is the
cheapest "DAG-free stays DAG-free" tripwire. The `hasattr` gate is the
illegal-states-unrepresentable structural assertion (the user's explicit ask 3a).

### (b) The per-mesh-shape DAG cache produces byte-identical graphs (no behavior change)

The cache move must NOT change a single graph. PROOF by inheritance: the
`window≡full` SWEEP + MATVEC oracles + the sweep-graph `window≡full` d1/d2/d3 all
stay byte-identical (`array_equal`) AFTER the move, because the DAG is the SAME
`SweepDependencyGraph.from_cartesian(spatial_shape)` — only its OWNER changes
(mesh → family accessor). Add ONE explicit cache-correctness assertion:

```python
@pytest.mark.foundation
def test_dag_cache_is_keyed_on_mesh_shape_and_byte_identical():
    """[L0 structural] Two meshes of the SAME spatial shape share byte-identical DAGs;
    the cached graph equals a freshly-built one (the move is a relocation, not a rebuild)."""
    repA = MovingFrontierWindow(_build_mesh("cart2d", nx=8, ny=6))
    repB = MovingFrontierWindow(_build_mesh("cart2d", nx=8, ny=6))   # same shape
    for label in repA._sweep_graphs:                 # the family's accessor
        gA = repA._sweep_graphs[label]
        gB = repB._sweep_graphs[label]
        # The DAG level structure + cell-index arrays are byte-identical.
        np.testing.assert_array_equal(gA.visit_order, gB.visit_order)
    # And the cache hit equals a fresh from_cartesian build (no behavior drift).
    fresh = SweepDependencyGraph.from_cartesian((8, 6))
    for label in repA._sweep_graphs:
        np.testing.assert_array_equal(repA._sweep_graphs[label].visit_order,
                                      fresh[label].visit_order)
```

(Adapt `visit_order` / the level-array accessor to the graph's actual public
attribute — read `sweep_graph.py` at implementation; the load-bearing assertion is
`array_equal` on the per-octant level/cell-index structure.) The PRIMARY proof
remains the window≡full oracles staying green — this test is the explicit
DAG-layer cross-check.

### (c) The `geometry.py:343/350` None-slot is GONE

```python
@pytest.mark.foundation
def test_no_sweep_graphs_none_slot_in_geometry():
    """[L0 structural] The curvilinear `self.sweep_graphs = None` slot is retired.

    geometry.py:343 (cylindrical) + :350 (spherical) set the None-slot today —
    the illegal-state smell. S6.4 removes them (the DAG is not a mesh concern)."""
    import inspect
    from orpheus.sn import geometry
    src = inspect.getsource(geometry)
    if "sweep_graphs" in src:
        pytest.fail(
            "geometry.py still mentions `sweep_graphs` — the curvilinear None-slot "
            "(was :343/:350) must be GONE; the DAG is owned by _DAGWavefront."
        )
```

(This is the cleanest retirement-audit assertion for the user's ask 3c. It also
catches the Cartesian `geometry.py:1377` build site — which must ALSO move to the
family accessor. If any legitimate `sweep_graphs` mention survives in geometry.py at
S6.4, this fails and the move is incomplete.)

---

## 4. The "one walk" discriminating test (the structural payoff)

Analogous to the S6.5 one-instance test ([[s6-relayering-verification]] §4): assert
`sweep` and `loss_action` go through the SAME `_OctantWalk2D` code, not two copies.
It MUST FAIL on current 3×-duplicated code and PASS after the carve — a genuine
tripwire, NOT a tautology. Two designs; ship the SPY (it directly observes the
"one walk" claim).

**File (drop-in, exact path):** `tests/sn/operators/test_one_octant_walk.py`

```python
r"""S6.4 — the 'one walk' discriminating test (#222).

The structural payoff of the unified kernel-parameterized octant walk: a SWEEP
((L+C)⁻¹q) and a LOSS_ACTION ((L+C)ψ, the matvec) on the SAME 2-D Cartesian mesh
MUST exercise the SAME ``_OctantWalk2D`` traversal code — they fork only at the
cell kernel (solve vs apply) + the emit policy, NOT in a duplicated octant frame.

NEGATIVE PRE-CONDITION (why this is a tripwire, not a tautology): pre-S6.4 the
octant orchestration is written THREE times in lockstep — ``sweep_octant_group``
(sweep), ``MovingFrontierWindow.loss_action`` (matvec), ``ScanMarch.loss_action``
(matvec) — so there is NO shared ``_OctantWalk2D`` for both directions to hit. The
test therefore FAILS at pre-S6.4 HEAD (the symbol does not exist / the two
directions touch different code objects) and PASSES after S6.4.

``-O``-safe: uses ``pytest.fail`` (fires under ``python -O``), NEVER a bare assert
(vv Mode 8).
"""
import numpy as np
import pytest

from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import CollisionOperator, StreamingOperator
from tests.sn._fixtures... import build_2d_cartesian_loss_operands  # reuse the affine-golden builder


@pytest.mark.foundation
def test_sweep_and_loss_action_hit_one_octant_walk(monkeypatch):
    """[L0 structural] sweep and loss_action exercise the SAME _OctantWalk2D._interior_walk."""
    # Import the shared walk; pre-S6.4 this import FAILS (symbol absent) → the test
    # fails loudly, recording the gap. Post-S6.4 it resolves.
    try:
        from orpheus.sn.loss_representation import _OctantWalk2D
    except ImportError:
        pytest.fail(
            "S6.4 NOT satisfied: `_OctantWalk2D` (the shared sweep+matvec octant "
            "walk) does not exist — the octant frame is still duplicated 3× in "
            "lockstep (sweep_octant_group + MovingFrontierWindow.loss_action + "
            "ScanMarch.loss_action). The one-walk unification (#222 S6.4) is open."
        )

    hits: list[str] = []
    real_walk = _OctantWalk2D._interior_walk
    def spy(self, *a, **k):
        hits.append("walked")
        return real_walk(self, *a, **k)
    monkeypatch.setattr(_OctantWalk2D, "_interior_walk", spy)

    L, C = build_2d_cartesian_loss_operands(ng=2, lvl=4)   # het-aniso ≥2G, the discriminating config
    A = L + C
    psi = _random_composite(A, np.random.default_rng(2026))

    # (1) the matvec direction goes through the shared walk:
    hits.clear()
    _ = L.apply(psi)                       # → loss_action → _OctantWalk2D._interior_walk
    if not hits:
        pytest.fail(
            "loss_action did NOT route through _OctantWalk2D._interior_walk — the "
            "matvec still uses a private duplicated octant frame (one-walk open)."
        )

    # (2) the sweep direction goes through the SAME shared walk:
    hits.clear()
    _ = A.solve(_random_source(A))         # → sweep → _OctantWalk2D._interior_walk
    if not hits:
        pytest.fail(
            "sweep did NOT route through _OctantWalk2D._interior_walk — the sweep "
            "still uses sweep_octant_group's private octant frame (one-walk open)."
        )
    # Both directions hit the ONE walk → L21 (matvec ≡ sweep) is a CODE FACT.
```

**Construction notes for the implementer:**
- The SPY observes the CLAIM directly: "do sweep AND matvec both reach
  `_OctantWalk2D._interior_walk`?" A tautology would spy on a method both already
  call by construction; here the LHS comes through `L.apply` (matvec door) and the
  RHS through `A.solve` (sweep door) — exactly the "are the two frames one frame"
  question. Pre-S6.4 the answer is structurally NO (3 frames, verified §0).
- **2-D config is mandatory** — the duplication is the 2-D octant frame
  (`sweep_octant_group` vs the matvec `loss_action`s). Use the SAME het-aniso ≥2G
  `build_2d_cartesian` operands as the affine golden so the walk is exercised on a
  non-flat, non-degenerate case (a 1G flat box makes the pure-z branch / interior
  recurrence trivial — vv §H1/§H2).
- **Stage as `xfail`→xpass:** at sub-step (a) (before `_OctantWalk2D` exists), land
  it `@pytest.mark.xfail(strict=True, reason="S6.4 extracts _OctantWalk2D shared by
  sweep+matvec; not yet cut")` keyed on the `ImportError`/no-hit `pytest.fail`. At
  sub-step (b) (when `sweep_octant_group` comes into the walk and both directions
  route through it), the `xfail` FLIPS to xpass — remove the marker IN the (b)
  commit. The recorded-gap discipline ([[scan-march-verification]] G3.b).
- **A simpler STRUCTURAL alternative** (if a spy is fragile against the
  cached/lazy frontier-plan construction): assert `sweep_octant_group` is GONE and
  both `MovingFrontierWindow.sweep` + `.loss_action` reference `_OctantWalk2D`:
  `inspect.getsource(...)` contains `_OctantWalk2D` for BOTH, AND
  `not hasattr(orpheus.sn.sweep, "sweep_octant_group")`. This is `-O`-safe and
  retirement-audit-aligned but observes the carve INDIRECTLY (source mention) rather
  than the runtime "one walk" fact. **Ship the SPY as primary; the source-structural
  as a cheap retirement-audit companion.**

---

## 5. The shared-walk OUTPUT-identity pin (the cheap relocation tripwire per sub-step)

Beyond the end-to-end oracles, add a CHEAP per-sub-step pin that the body did NOT
change when it moved into `_OctantWalk2D` — a direct equivalence at the
representation boundary on a het-aniso ≥2G config (Mode-9: NOT the degenerate box).

**File:** extend `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py` (the
oracle home) OR a new `tests/sn/operators/test_octant_walk_relocation.py`.

```python
@pytest.mark.foundation
@pytest.mark.parametrize("nx,ny,lvl,ng,bc", [(12, 7, 6, 2, "reflective"), (5, 9, 4, 4, "reflective")])
def test_window_matvec_byte_identical_across_octant_walk_extraction(nx, ny, lvl, ng, bc):
    """[L0] MovingFrontierWindow.loss_action output is byte-identical to the
    FullFieldWavefront oracle AFTER the octant frame relocates into _OctantWalk2D.

    The relocation is a pure code move (octant projection + boundary block →
    shared walk); the matvec output MUST stay array_equal. NON-SQUARE catches an
    x↔y swap the shared walk could introduce."""
    L, psi = _build_2d_loss_and_probe(nx, ny, lvl, ng, bc)   # het-aniso ≥2G
    win  = MovingFrontierWindow(L.sn_mesh).loss_action(L, psi)
    full = FullFieldWavefront(L.sn_mesh).loss_action(L, psi)
    np.testing.assert_array_equal(win.bulk.values, full.bulk.values, err_msg="bulk (L+C)ψ")
    for face in L.sn_mesh.trace.face_names:
        np.testing.assert_array_equal(
            win.boundary.face_view(face), full.boundary.face_view(face),
            err_msg=f"boundary defect {face}",
        )
```

This is the relocation-identity pin: the move is byte-preserving, so `win ≡ full`
bit-identical at every sub-step that touches the octant frame. It is the
representation-boundary twin of the end-to-end `test_matvec_window_equals_full_field`
(same assertion, cheaper to run per sub-step). (If the existing
`test_matvec_window_equals_full_field_end_to_end` already drives exactly these two
representations through `op.apply`, this can be subsumed by it — confirm the existing
oracle hits `MovingFrontierWindow` vs `FullFieldWavefront` and not the operator's
default; if it does, cite it instead of adding this.)

---

## 6. Anti-degradation gate — the kernel parameter MUST stay a KERNEL OBJECT

CONSTRAINT (the brief; the boolean-flag anti-pattern, `coding-elegance` Smell #3):
`_OctantWalk2D` is parameterized by a CELL-KERNEL OBJECT (the cell operation
`apply_windowed`/`residual_windowed`/the row-march scan) + an EMIT policy, **NOT a
boolean `is_solve` flag**. Gate the design so the carve cannot degrade into a flag.

**File:** `tests/sn/operators/test_one_octant_walk.py` (same file as §4).

```python
@pytest.mark.foundation
def test_octant_walk_is_kernel_parameterized_not_boolean():
    """[L0 structural] _OctantWalk2D forks on a KERNEL OBJECT + EMIT policy,
    NOT a boolean is_solve/is_apply flag (coding-elegance Smell #3).

    A boolean dispatch (`if solve: cell_kernel_batch else residual_kernel_batch`)
    inside the walk is the anti-pattern: it re-introduces the twin-path branch the
    carve eliminates. The walk's interior MUST receive its cell operation as a
    callable/object parameter."""
    import inspect
    from orpheus.sn.loss_representation import _OctantWalk2D
    src = inspect.getsource(_OctantWalk2D)
    # No boolean is_solve/is_apply branch inside the shared walk.
    for smell in ("is_solve", "is_apply", "is_matvec"):
        if smell in src:
            pytest.fail(
                f"_OctantWalk2D contains a boolean `{smell}` dispatch — the carve "
                "degraded into the boolean-flag anti-pattern (Smell #3). The cell "
                "operation MUST be a kernel OBJECT/callable parameter, not a flag."
            )
```

This is a STRUCTURAL guard against the specific degradation the brief names. It is
necessarily a source-inspection test (the anti-pattern is a code shape, not a
runtime value) — `-O`-safe (`pytest.fail` + `inspect.getsource`). Pair it with the
elegance-enforcer review (§7) which is the authoritative judge of "kernel object vs
flag"; this test is the mechanical tripwire that the review's verdict does not
silently regress.

---

## 7. Staging — each sub-step independently bit-identity-gated

RECOMMENDED ORDER (cheapest-first / deepest-unification-progressive; each
independently committable + gated):

| Sub-step | What moves | Gates that MUST stay green | NEW test(s) | bit-id vs principled | `-O` note |
| -------- | ---------- | -------------------------- | ----------- | -------------------- | --------- |
| **(a) Extract `_OctantWalk2D` for the TWO MATVEC variants** (window + scanmarch): the shared octant projection + pure-z + face strings + boundary block; `MovingFrontierWindow.loss_action` + `ScanMarch.loss_action` delegate, providing only their `_interior_walk` (residual_windowed | row-march scan). | the matvec octant frame + the O.4b boundary block (×2 → ×1) | **window≡full MATVEC oracle** (the output-identity proof); ScanMarch G2.c (nulp); affine-carve golden; `loss_action` convention pin (`cart2d_2g`) | §5 relocation-identity pin; §6 anti-boolean guard | **A2D-1 RETIRES** (§2) — or retargets to `_OctantWalk2D._interior_walk` + regen in-commit; everything else **bit-identity** (relocation byte-preserving) | §5/§6 `-O`-safe; retire A2D-1's `getsource` ref |
| **(b) Bring `sweep_octant_group` + `_sweep_2d_scanmarch` into `_OctantWalk2D`** (the deepest L21 unification): the SWEEP octant frame routes through the SAME walk, forking at `cell_kernel_batch` (solve) + the ψ/scalar/moment emit. | the sweep octant frame (`sweep.py:796`, `:1276`) into the shared walk | **window≡full SWEEP oracle** (the output-identity proof for the sweep direction); ScanMarch G2.c sweep; affine-carve golden; SI≡Krylov≡k_inf | **§4 one-walk SPY FLIPS `xfail`→xpass** (the structural payoff) | **bit-identity** (relocation byte-preserving) | §4 `pytest.fail` (`-O`-safe); the `xfail` removal is the visible green flip; retire the `sweep.py:1328` IOU marker |
| **(c) DAG-ownership move** (mesh → `_DAGWavefront` family, cached per mesh-shape): `geometry.py:343/350/1377` `sweep_graphs` → the family accessor; `loss_representation.py:524/672` + `sweep.py:910/1533` re-point. | the DAG off the mesh | **window≡full SWEEP + MATVEC** (the cache is byte-identical); sweep-graph window≡full d1/d2/d3; affine-carve golden; SI≡Krylov≡k_inf | **§3 (a) mesh-no-sweep_graphs + (b) cache byte-identical + (c) geometry None-slot gone** | **bit-identity** (relocation; same `from_cartesian` graph, new owner) | §3 tests `-O`-safe; retirement audit: grep `sweep_graphs` in geometry.py → expect zero |
| **(d) Fold `FullFieldWavefront`** into the shared boundary block (the THIRD octant frame + boundary block, `loss_representation.py:611`): the d-generic interior walk stays the full-cochain `_interior_walk`, the octant projection + boundary block come from `_OctantWalk2D`. | the oracle's octant frame + boundary block (×3 → ×1 total) | **window≡full SWEEP + MATVEC** (the oracle MUST stay byte-identical — it is the reference); sweep-graph window≡full d3 (oracle d-generic); affine-carve golden | (none new — consolidation) | **bit-identity** (the oracle output is the reference; if it drifts, the carve is wrong) | retirement audit: the THREE lockstep IOU markers (`loss_representation.py:479/:808` + `sweep.py:1328`) GONE; ONE `_OctantWalk2D` frame; one boundary block |

**Why this order:** (a) is cheapest (two matvec frames already in one module,
side-by-side, docstring-flagged "IDENTICAL"), fires the A2D-1 decision early, and is
gated by the MATVEC output oracle alone. (b) is the deep unification (crosses the
sweep.py↔loss_representation.py boundary) and is where the one-walk SPY flips — do it
on top of a clean (a). (c) is orthogonal structure (mesh → family) — independently
gated by the cache byte-identity, can interleave after (a) or (b). (d) folds the
oracle last (it is the REFERENCE — fold it only when the production frames are
already shared, so the oracle's byte-identity is the final consolidation proof, not a
moving target during the carve).

**Retirement audit (elegance-enforcer deliverable, gated at (d)):** ZERO "edit both
in lockstep" markers survive (grep the three strings → expect zero); ONE
`_OctantWalk2D` octant frame; ONE O.4b boundary block (was ×4: 3 matvec + the sweep
shed); `sweep_octant_group` + `_sweep_2d_scanmarch`'s private frames GONE; A2D-1
retired (or retargeted); `geometry.py` free of `sweep_graphs`. The §4 one-walk SPY
green is the construction-level proof the lockstep duplication is closed.

---

## 8. Mode-8 / Mode-9 discipline (explicit)

- **Mode-8 (`-O` strips bare `assert`).** EVERY new S6.4 tripwire is `-O`-safe:
  §3/§4/§5/§6 use `pytest.fail` + `np.testing.assert_array_equal` (function calls
  that fire under `-O`). Grep `^\s*assert ` in every new S6.4 test file; expect
  zero. The neighbourhood was already cleaned at S6.3a (`test_g_adjoint_reciprocity`
  migrated to `pytest.fail`, `ee91db4`) — confirmed this session (`:210/:222/:248/
  :286/:303` all `pytest.fail`). The collection advisory ("assert statements are
  not executed … python -O") fires because OTHER files in the broad collection still
  carry bare asserts — NOT the S6.4 gate set. Do NOT introduce a bare assert.
- **Mode-9 (FP-invariance only in a degenerate regime).** S6.4 is a RE-LAYERING, not
  a splitting/accelerator — the converged FP does NOT change (forks no default). The
  classic Mode-9 splitting trap does not directly apply. BUT every equivalence /
  relocation-identity gate MUST run on an ANISOTROPIC + HETEROGENEOUS + ≥2G config
  (NON-SQUARE 12×7/5×9, LS-4, 2G het aniso — already the oracle + G2.c config), NOT
  the isotropic-reflective box: the octant interior recurrence + the pure-z branch +
  the x↔y face mapping ALL degenerate on a 1G flat square box (vv §H1/§H2 + the
  [[phase5a_5c_angular_windowing_lessons]] SQUARE-snapshot x↔y blindness). The
  NON-SQUARE shapes in the anchor set are LOAD-BEARING — they catch the
  variable-swap (mu_x↔mu_y) the shared walk could introduce (vv failure mode #2). If
  S6.4 ever flips a default (it does not, by design), the FP-invariance gate inherits
  the full Mode-9 discipline (diagonal cubature + vacuum/streaming) from
  [[scan-march-verification]] §G4.

---

## 9. Deselect / standing reds (inherited)
- **DESELECT** `tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff`
  (#212 `continuous_get` hang).
- Held reds untouched by the 2-D-Cartesian octant-frame carve: **#206** cyl-matvec
  (deferred curvilinear adjoint — S6.4 does NOT touch the curvilinear 1-D
  `CumprodScan` path; stays `xfail`); **#195** MMS@160.
- The 2-D Cartesian ADJOINT deferral (`_DAGWavefront.loss_action_transpose` raises
  `NotImplementedError`, `loss_representation.py:411-425`) MUST be PRESERVED through
  S6.4 — the shared `_OctantWalk2D` is FORWARD-only; the multi-D adjoint stays the
  deferred raise (NEVER a silent wrong answer). If a `loss_action_transpose` shared
  walk is added later it is a SEPARATE carve with its own gate.

## 10. Self-improvement triggers (fired BEFORE delivery)
- **NEW failure mode?** S6.4 introduces NO new failure mode beyond the
  `vv-principles` table. The "three octant frames in lockstep" defect is
  `coding-elegance` Smell #16 / Pattern-2 duplication (already cataloged); the
  boolean-flag degradation risk is Smell #3 (gated by §6). The x↔y SQUARE-snapshot
  blindness is the [[phase5a_5c_angular_windowing_lessons]] TRAP (already a memory
  note). NO skill-table append required.
- **Plan rejection?** N/A (pre-implementation). If qa/user rejects, log the
  counter-example per the AGENT.md self-improvement directive in `feedback_*.md`.

## 11. Pre-reads cross-check (file:line, VERIFIED this session @ 7bbeb5f)
- `orpheus/sn/loss_representation.py:451` — `MovingFrontierWindow.loss_action` (matvec octant frame #1; `graph.residual_windowed` `:534`; O.4b boundary block `:547-558`).
- `orpheus/sn/loss_representation.py:611` — `FullFieldWavefront.loss_action` (matvec octant frame #2, d-generic; `graph.residual` `:675`; O.4b block `:689-699`).
- `orpheus/sn/loss_representation.py:779` — `ScanMarch.loss_action` (matvec octant frame #3, row-march; `_x_scan_faces` `:877`; O.4b block `:897-909`, docstring "IDENTICAL to MovingFrontierWindow.loss_action").
- `orpheus/sn/loss_representation.py:479` + `:808` — the "edit both in lockstep" IOU markers (matvec, S6.4 retires).
- `orpheus/sn/loss_representation.py:389-425` — `_DAGWavefront` base + the deferred 2-D adjoint raise (`:411-425`, preserve).
- `orpheus/sn/sweep.py:796` — `sweep_octant_group` (sweep octant frame; `graph.apply_windowed` `:936`; outflow shed `:956-957`).
- `orpheus/sn/sweep.py:1276` — `_sweep_2d_scanmarch` (sweep row-march; docstring `:1323-1335` "edit both in lockstep" with `sweep_octant_group`, IOU at `:1328`).
- `orpheus/sn/sweep.py:1451` — `_sweep_full_field` (oracle sweep counterpart).
- `orpheus/sn/spatial/diamond.py:273` — `cell_kernel_batch` (SOLVE kernel); `:331` — `residual_kernel_batch` (APPLY kernel). The fork point.
- `orpheus/sn/geometry.py:343` (cyl) + `:350` (sph) — `self.sweep_graphs = None` (the None-slot S6.4 removes); `:1377` — the Cartesian DAG build (moves to the family accessor).
- `tests/sn/operators/test_streaming_operator.py:1321/:1333` — A2D-1 `EXPECTED_SHA256` on `getsource(MovingFrontierWindow.loss_action)` (§2 retire/retarget target).
- `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py:77/:101` — window≡full SWEEP + MATVEC oracles (4 cases incl. NON-SQUARE 12×7/5×9).
- `tests/sn/sweep/cartesian_2d/test_scan_march_equivalence.py:93/:118/:168` — ScanMarch G2.c (sweep/moment/residual ≡ oracle, nulp).
- `tests/sn/sweep/core/test_sweep_graph_window_equivalence.py:197/:219` — sweep-graph window≡full solve/residual d1/d2/d3 (synthetic admission).
- `tests/sn/solve/test_affine_carve_bit_identity.py:148` — the end-to-end sha256 golden (3 params).
- `tests/sn/operators/test_loss_action_convention.py` — `[slab_2g]`/`[cart2d_2g]` the `(L+C)ψ` + `−C` glue pin (§1).
- `tests/sn/operators/test_g_adjoint_reciprocity.py:210+` — the curvilinear reciprocity sentinel (Mode-8-fixed at S6.3a; stays as regression, NOT re-carved by S6.4).

---

# ADDENDUM (2026-06-11): sub-steps (e) + (f)

**Trigger:** user-approved S6.4 AMENDMENT (LOCKED 2026-06-11; the native-place
architecture audit), recorded in `.claude/plans/sn_sweep_strategy.md` §"S6.4 —
AMENDMENT" (`:719-763` of the worktree plan). Worktree HEAD advanced `7bbeb5f`→
`9cfb7a1` since the body above. Two NEW sub-steps were appended to the staged
(a)–(d). Sub-steps (a)–(d) above stay EXACTLY as planned — with the ONE naming
correction the amendment locks: the shared walk is **`_OctantWalk` (d-GENERIC from
birth, per-axis sign/face/inflow tuples over `mesh.ndim`), NOT `_OctantWalk2D`**
(folding the d-generic full-field oracle into a 2-D-hardcoded walk at (d) would
regress C3). Everywhere the body above writes `_OctantWalk2D`, read `_OctantWalk`;
its `_interior_walk` is byte-identical at d=2, so every gate above transfers
unchanged (the SPY §4, the anti-boolean guard §6, the relocation pin §5).

The audit's core finding: the solve/apply direction fork bottoms out in
`orpheus/sn/spatial/diamond.py` — **VERIFIED this session**: `update_batch`
(`:453`, SOLVE) / `residual_batch` (`:513`, APPLY), the pure, storage-free,
d-generic cell kernels (the docstrings `:370`/`:540` already flag the gather +
`denom` as "IDENTICAL between update_batch (solve) and residual_batch (apply)").
That fork is then RE-SPELLED at three altitudes: the octant frame (the (a)–(d)
target), **the graph level-walks** (sub-step (e)), and the cell storage adapter
(thin already). The module geography (bodies living apart from their owner +
the load-time import cycle) is sub-step (f).

## E.0 — Pre-reads VERIFIED @ HEAD `9cfb7a1` (the (e)/(f) ground truth)

| Symbol | file:line | Role |
| ------ | --------- | ---- |
| `SweepDependencyGraph.apply` | `sweep_graph.py:682` | full-field SOLVE level-walk (→ `update_batch` `:755`) |
| `SweepDependencyGraph.residual` | `sweep_graph.py:770` | full-field APPLY level-walk (→ `residual_batch` `:825`) |
| `SweepDependencyGraph._make_slice` | `sweep_graph.py:644` | the SHARED level/cell-index packet builder for apply+residual (`:657` "Shared by apply…and residual") |
| `SweepDependencyGraph.apply_windowed` | `sweep_graph.py:849` | windowed SOLVE level-walk (rebuilds the `_MovingFrontier` seed/incoming/emit/shed loop) |
| `SweepDependencyGraph.residual_windowed` | `sweep_graph.py:929` | windowed APPLY level-walk (rebuilds the SAME frontier loop) |
| `SweepDependencyGraph.from_cartesian` | `sweep_graph.py:445` | DAG construction (d-generic); **`for_shape` does NOT exist yet** |
| `CellUpdate` Protocol | `spatial/cell_update.py:467` | the override point (`update`/`residual` per-cell `:505`/`:553`) |
| `CellUpdateBase.update_batch` / `.residual_batch` | `spatial/cell_update.py:714`/`:751` | the batched kernel pair (DD overrides in `diamond.py:453`/`:513`) |
| `_x_scan_faces` / `_scanmarch_row` | `sweep.py:1186`/`:1250` | the scan-family bodies → (f) relocate to `spatial/scan.py` |
| `_sweep_2d_scheduled` | `sweep.py:960` | the schedule driver → (f) relocate to `loss_representation.py` |
| `_sweep_1d_unified` / `transport_sweep` | `sweep.py:315`/`:104` | → (f) relocate to `loss_representation.py` |
| `sweep_graphs` (mesh-owned today) | `geometry.py:1377` (Cartesian build) + `:343`/`:350` (curvilinear `=None`) | (c) moves off mesh; (f) re-homes the cache via `for_shape` |

**CLAIM VERIFIED (the user's ask for (e)): the runtime consumers of
`update_batch`/`residual_batch` are EXACTLY `sweep_graph.py:755` (apply) and `:825`
(residual).** No other PRODUCTION call site exists — the only other `update_batch`/
`residual_batch` mentions in `orpheus/` are docstrings (`sweep.py:1107`,
`loss_representation.py:501`, `diamond.py` def-site docstrings, and the
**`cell_update.py:777-779` round-trip-contract example is a DOCSTRING, NOT a
runtime consumer** — confirmed by reading the surrounding `residual_batch`
docstring). TEST call sites: `tests/sn/sweep/core/test_cell_update_batch.py` (15+
direct `DiamondDifference().update_batch(...)` / `.residual_batch(...)` calls),
`test_sweep_graph.py`, `test_sweep_graph_nd_admission.py`,
`test_streaming_operator.py`. So the brief's "(e) makes the CellUpdate override
point the pure kernel pair instead of update_batch/residual_batch" is sound — the
batch methods are consumed only by graph.apply/residual in production, and the
direct-call tests REWIRE to the new pure kernel pair (retirement = test migration).

## E.1 — Sub-step (e): collapse the graph's 4 level-walks → 2 direction-parameterized walks

**What moves.** Today `sweep_graph.py` carries FOUR level-walks:
`apply`(:682)/`residual`(:770) (the full-field pair, already sharing `_make_slice`)
and `apply_windowed`(:849)/`residual_windowed`(:929) (the window pair, EACH
rebuilding the `_MovingFrontier` seed/incoming/emit/shed level loop). After (e):
**TWO buffer-policy walks** — `full` (the `_make_slice` cochain loop) and
`windowed` (the `_MovingFrontier` loop) — **each parameterized by a DIRECTION
object** = (cell kernel + per-level emit policy). The direction objects are the
SAME emit-policy objects designed at sub-step (a) (the solve emit = ψ angular/
moment write; the apply emit = `LpC` bulk + O.4b boundary defect) — **reuse, NOT
re-mint** (Pattern 2: single source of truth; the (a) emit policy is the producer).
Also: the `CellUpdate` Protocol override point for future Step/LD closures becomes
the **pure kernel pair** (`update_batch`/`residual_batch` → the cell algebra a
strategy supplies), with storage/gather/scatter handled ONCE in the walk above.

### E.1.1 — Per-sub-step gate table (same format as §7)

| What moves | Anchors that MUST stay green | NEW test(s) | bit-id class | `-O` note |
| ---------- | ---------------------------- | ----------- | ------------ | --------- |
| The 4 level-walks (`apply`/`residual`/`apply_windowed`/`residual_windowed`) → 2 direction-parameterized walks (`full`,`windowed`); the (a) emit-policy objects threaded as the direction param; `CellUpdate` override point → the pure kernel pair | **window≡full SWEEP + MATVEC** (`test_2d_full_field_oracle.py`, `array_equal`, NON-SQUARE 12×7/5×9 — the BOTH-directions proof the 2 collapsed walks reproduce solve AND apply); **sweep-graph window≡full d1/d2/d3** (`test_sweep_graph_window_equivalence.py` — the GRAPH-LAYER proof the 4→2 collapse is byte-identical at the layer it touches, incl. synthetic d=3); **ScanMarch G2.c** (nulp — the scan path threads the same emit); **affine-carve golden** sha256; **SI≡Krylov≡k_inf** | §E.1.2 emit-order pin (per-level einsum accumulation byte-identity, d=2 + d=3); §E.1.3 left-fold-order guard (kernel left-fold UNTOUCHED); §E.1.4 CellUpdate-override-is-kernel-pair structural test | **bit-identity** (relocation byte-preserving; the 4→2 collapse is a re-org of the SAME level loop + SAME kernel calls, no FP-reduction-tree change) | all `pytest.fail`/`assert_array_equal` (Mode-8-safe); grep `^\s*assert ` in new files → zero |

### E.1.2 — Is a NEW discriminating/structural test warranted for (e)? YES — but a TIGHT one

The output oracles (window≡full both directions + sweep-graph window≡full d1/d2/d3)
**suffice for the value/behavior claim** — they are `array_equal`, and the 4→2
collapse changes no value. BUT three things could SILENTLY break in a way the
end-to-end oracles are slower to localize, so add cheap GRAPH-LAYER pins:

1. **Per-level emit REORDER risk (THE real (e) hazard).** The scalar/moment emit
   `+=` accumulates across anti-diagonal levels (`moment_buf[...] += einsum(...)`)
   AND the solve emit writes ψ per level. **This accumulation order is
   bit-id-load-bearing WITHIN a schedule** — if the 4→2 collapse re-orders the
   per-level emit (e.g. iterates levels in a different sequence, or batches the
   emit differently), the moment `+=` reduction tree changes → `array_equal`
   FAILS (correctly), but the failure surfaces as a whole-suite red far from the
   cause. PIN it at the graph layer: a foundation test that the `full` and
   `windowed` walks emit byte-identical moment buffers on a het-aniso ≥2G NON-SQUARE
   config, AND that the level ITERATION ORDER (`visit_order` per octant) is
   unchanged from a fresh `from_cartesian` build (so the collapse re-uses, never
   re-derives, the level structure). This is the moment-mode analogue of the
   [[phase5a_5c_angular_windowing_lessons]] THE-TRAP catcher — a scalar-only check
   is blind to ℓ/m drift, so the pin MUST assert the FULL moment tensor (incl.
   ℓ≥1), not just `moments[0,0]`.

   ```python
   @pytest.mark.foundation
   @pytest.mark.parametrize("shape", [(12, 7), (5, 9), (3, 2, 3)])  # NON-SQUARE + synthetic d=3
   def test_collapsed_walks_emit_byte_identical_moments_and_share_visit_order(shape):
       """[L0] The 4→2 level-walk collapse preserves the per-level emit accumulation
       order (moment += is FP-order-load-bearing within a schedule) AND re-uses the
       level structure (never re-derives it)."""
       g = SweepDependencyGraph.from_cartesian(shape)
       # full vs windowed must agree byte-for-byte on the FULL moment tensor (ℓ≥1):
       m_full = _run_full_walk_moments(g, _het_aniso_2g_probe(shape))
       m_win  = _run_windowed_walk_moments(g, _het_aniso_2g_probe(shape))
       np.testing.assert_array_equal(m_full, m_win, err_msg="moment emit reorder")
       # level iteration order is the SAME object as a fresh build (re-use, not re-derive):
       fresh = SweepDependencyGraph.from_cartesian(shape)
       for label in g.octant_labels:
           np.testing.assert_array_equal(g.visit_order(label), fresh.visit_order(label))
   ```
   (Adapt `visit_order`/`octant_labels`/the moment-run helpers to the actual graph
   surface — read `sweep_graph.py` at impl; the load-bearing assertions are
   `array_equal` on the FULL moment tensor + the per-octant level structure.)

2. **Left-fold order in the KERNELS must NOT be touched.** The cell kernels
   (`update_batch`/`residual_batch`, the `α`-recurrence + the WDD `out=2ψ̄−in` +
   the per-Legendre `R·Λ·M` left-fold) are the DD reduction tree the affine golden
   inherits from. (e) collapses the LEVEL-WALKS, NOT the kernels — the kernel
   bodies are pure and storage-free and STAY in `diamond.py`. Guard that (e) does
   not "helpfully" inline or re-order the kernel fold while parameterizing the
   walks: a foundation source-structural pin that `diamond.py::update_batch` /
   `residual_batch` are unchanged (an `inspect.getsource` sha256 on the TWO kernel
   bodies — the legitimate "this body is the FP-reduction tree of record" pin,
   the genuine STRUCTURAL-reference exception to aggressive-retirement, NOT the
   A2D-1 same-math-text case). If the implementer touches the kernel fold, the
   sha256 trips AND the affine golden trips — paired loud failure. (If a kernel
   hash is judged too brittle, the affine golden + window≡full both-directions
   alone catch any fold change as a value drift; the hash is the cheap localizer.)

3. **The moment-mode `+=` principled-equivalence boundary.** Moment output is a
   principled-equivalence (nulp) carve where it differs from post-projection
   (the [[s6-relayering-verification]] fuller-view oracle) — but WITHIN the (e)
   collapse the moment emit is bit-identical (same einsum, same level order). So
   (e)'s moment pin is `array_equal` (NOT nulp): the collapse must not introduce
   ANY moment FP drift. The nulp boundary stays where it already lives (the
   `solve_moments ≡ solve+project` fuller-view oracle), untouched by (e).

### E.1.5 — What the (e) output oracles DO suffice for (no new test needed)

The "4 old graph methods GONE / 2 new walks shared by both directions" CLAIM does
NOT need its own SPY beyond the §4 one-walk SPY — the §4 SPY already observes that
sweep AND matvec route through the shared `_OctantWalk._interior_walk`, and the
walk's level loop IS one of the 2 collapsed walks. ADD only a cheap retirement-audit
structural assertion (it belongs to the (e) retirement, not a behavior gate):
`apply_windowed`/`residual_windowed` are GONE as separate methods (the windowed walk
is now ONE method parameterized by direction) — `not hasattr(SweepDependencyGraph,
"apply_windowed")` and `not hasattr(..., "residual_windowed")` (or, if the names are
kept as thin direction-bound wrappers, assert their bodies are ≤2 lines delegating to
the parameterized walk via `inspect.getsource`). `-O`-safe (`pytest.fail`).

## E.2 — Sub-step (f): module geography — `sweep.py` DISSOLVES

**Relocations (VERIFIED targets @ `9cfb7a1`):**
- `_x_scan_faces`(`sweep.py:1186`) / `_scanmarch_row`(`:1250`) → `spatial/scan.py`
  (the scan family's "graph+walk together"; `spatial/scan.py` ALREADY hosts
  `ordinate_scan` and cites `_sweep_1d_unified` at `:41` — the natural home).
- the octant walk (`_OctantWalk`) + the kernel-parameterized `_sweep_2d_scheduled`
  schedule driver (`:960`) + `_sweep_1d_unified`(`:315`) + `transport_sweep`(`:104`)
  → `loss_representation.py`.
- DAG construction cache = `SweepDependencyGraph.for_shape(shape)` classmethod
  (cached per mesh-shape) accessed via the `_DAGWavefront` family accessor.
- `solver.py:344` re-points its `_sweep_2d_scheduled` import; `solver.py:55`
  re-points `transport_sweep`. `_reflect_outflow_into_inflow`(`solver.py:1428`) is
  ALREADY native in `solver.py` (a boundary concern) — stays.
- the **loss_representation↔sweep lazy-import cycle dissolves** (VERIFIED today:
  `loss_representation.py:135` imports the sweep bodies at module level;
  `sweep.py:239` lazily imports `default_for` back to break it — after (f) the
  bodies live IN `loss_representation.py`, no back-edge, no workaround).
- production consumers: `operator.py:2094` (`transport_sweep`), `solver.py:55/344`.
  ~25 test files import `from orpheus.sn.sweep` (VERIFIED grep count = 25, the
  brief's "~30" rounds this).

### E.2.1 — The (f) gate is PURE RELOCATION (no value claim; behavior frozen)

| Gate | Mechanism | Pass condition |
| ---- | --------- | -------------- |
| Clean collection under `-O` | `python -O -m pytest --collect-only tests/sn` | exit 0, zero ImportError/collection errors |
| **Zero `from orpheus.sn.sweep import` survivors** (grep-gate) | `grep -rn "from orpheus.sn.sweep import\|orpheus\.sn\.sweep\." orpheus/ tests/` (excl. `sweep.py` if a thin re-export shim is kept ≤1 merge cycle) | zero in `orpheus/` PRODUCTION; tests rewired to `spatial/scan` + `loss_representation` homes |
| Full anchor set green | window≡full SWEEP+MATVEC, sweep-graph window≡full d1/d2/d3, ScanMarch G2.c, affine-carve golden sha256, SI≡Krylov≡k_inf, loss_action convention pin, the §4 one-walk SPY (now imports `_OctantWalk` from its NEW home in `loss_representation.py`) | all green, byte-identical (relocation preserves FP-reduction tree exactly) |
| Sphinx `:func:`/`:mod:` roles updated | rebuild `sphinx-build docs docs/_build/html -W` | clean build, zero "py:func reference target not found" for relocated symbols |

### E.2.2 — The (f) RETIREMENT-AUDIT checklist (elegance-enforcer deliverable)

VERIFIED scope this session — the Sphinx surface is SUBSTANTIAL and MUST migrate
(the densest dangling-role risk):

- [ ] `sweep.py` DELETED (or reduced to a ≤1-merge-cycle thin re-export shim that
      the NEXT merge removes — per [[feedback_aggressive_retirement]] deprecation
      shims live one cycle only). The amendment's intent is full dissolution.
- [ ] `__all__` / module-level exports migrated to the new homes; no symbol exported
      from two modules.
- [ ] **Sphinx `:func:`/`:mod:` role migration (VERIFIED dangling targets):**
      `docs/theory/index_convention.rst` (~14 `transport_sweep` roles, `:444`–`:1347`),
      `docs/theory/boundary_conditions.rst` (`:1547` `:mod:orpheus.sn.sweep`, `:2030`
      `transport_sweep`, `:2031/:2055` `_sweep_1d_spherical`, `:2111` `_sweep_1d_cylindrical`),
      `docs/theory/operator_algebra.rst` (`:1512`/`:2109`/`:2128` `transport_sweep`+
      `_sweep_2d_wavefront`; `:5052`–`:5067` `_sweep_2d_scheduled`/`sweep_octant_group`/
      `_sweep_2d_wavefront`/`transport_sweep` thread-the-kernel block),
      `docs/theory/verification.rst` (`:46` `_sweep_1d_cumprod`),
      `docs/theory/structured_geometry.rst` (`:392` `orpheus/sn/sweep.py` file-path ref).
      Each `~orpheus.sn.sweep.X` → `~orpheus.sn.loss_representation.X` (walk/driver/
      `transport_sweep`/`_sweep_1d_unified`) or `~orpheus.sn.spatial.scan.X`
      (`_x_scan_faces`/`_scanmarch_row`/scan bodies). The `-W` Sphinx build is the
      mechanical gate: a stale `:func:` is a build error, NOT a silent dangle.
- [ ] In-`orpheus/` docstring cross-refs to `~orpheus.sn.sweep.*` updated (VERIFIED:
      `scattering.py:695`, `sweep_graph.py:7`, `operator.py:1087/1535/1808/2054`,
      `diamond.py:14/457`, `loss_representation.py:461/618/724/788`,
      `spatial/sweep_cache.py:76`, `spatial/scan.py:41`). These are `:func:` roles in
      docstrings — `-W` catches them too once Nexus/autodoc re-resolves.
- [ ] ~25 test files rewired to the new import homes (retirement = test migration,
      [[feedback_retirement_means_test_migration]]); no test imports a dead path.
- [ ] No `:file:orpheus/sn/sweep.py` path reference survives as a live pointer
      (`structured_geometry.rst:392`, `boundary_conditions.rst:2122` mention the
      file path in prose — update to the new module).
- [ ] The lazy-import cycle is GONE: `grep "from .loss_representation import"` in the
      new `loss_representation.py` shows NO back-edge to a sweep module (the cycle
      that motivated the lazy import no longer exists).

### E.2.3 — (f) needs NO new behavior test

(f) is byte-preserving relocation: the full anchor set (which already includes the
§4 SPY re-imported from the new home, the window≡full both-directions, the affine
golden) IS the behavior gate. The ONLY new artifacts are the grep-gate (zero
`from orpheus.sn.sweep import` in production) and the `-W` Sphinx build — both
mechanical, both `-O`-irrelevant (collection + build, not assertion). Do NOT add a
value test for a pure module move.

## E.3 — SEQUENCING CHECK: (e)/(f) AFTER (a)–(d), and the (c)↔(f) `for_shape` hazard

**(e) and (f) come AFTER (a)–(d)** — confirmed sound:
- (e) collapses the level-walks the (a) emit-policy objects parameterize; the (a)
  emit objects MUST exist first (Pattern 2 reuse). And (b) must have routed sweep
  into the shared walk before (e) generalizes the walk's level loop, else (e) has
  only the matvec half to collapse.
- (f) is the geography move; it relocates the `_OctantWalk` + driver that (a)–(d)
  built and (e) finalized. Moving a still-changing body is churn — (f) LAST.

**THE FLAGGED HAZARD — (c) and (f) BOTH touch the `sweep_graphs` / DAG-cache
symbol.** VERIFIED state @ `9cfb7a1`: **`SweepDependencyGraph.for_shape` does NOT
exist** (only `from_cartesian` at `:445`); the DAG is mesh-owned (`geometry.py:1377`
Cartesian build + `:343/:350` curvilinear `None`-slot). The original sub-step (c)
(§3 above) moves OWNERSHIP off the mesh into the `_DAGWavefront` family accessor
but does NOT, as written, mint the `for_shape` classmethod by name. The (f) bullet
introduces `SweepDependencyGraph.for_shape(shape)` as the cached construction the
family accessor calls. **If (c) moves the cache to an ad-hoc family-held dict and
(f) THEN refactors that into `for_shape`, the same DAG-cache surface is touched
twice — two relocations, two A2D-1-style regenerations of the §3(b) cache
byte-identity pin, double the churn the audit warns against.**

**RECOMMENDATION: (c) SHOULD already create `SweepDependencyGraph.for_shape` on the
graph class.** Move the construction (with the graph) AND the ownership (with the
family) in ONE step at (c): `for_shape(shape)` is the cached classmethod
(`@lru_cache`/dict keyed on `shape`, calling `from_cartesian`); the `_DAGWavefront`
family accessor calls `for_shape`. Then (f) only re-points imports and dissolves
the module — it does NOT touch the DAG-cache mechanism again. This collapses the
two touches into one, and the §3(b) `test_dag_cache_is_keyed_on_mesh_shape_and_byte_identical`
pin (adapt its accessor to `for_shape`) gates the SINGLE move. Update the §3
sub-step (c) gate to assert `for_shape` exists on the graph class AND the family
accessor routes through it (`hasattr(SweepDependencyGraph, "for_shape")` +
`not hasattr(SNMesh, "sweep_graphs")` for all coords). The brief's "DAG construction
cache = `SweepDependencyGraph.for_shape(shape)` … sub-step (c) unchanged,
better-homed" (plan `:756-757`) AGREES: home `for_shape` at (c), not (f).

## E.4 — Mode-8 / Mode-9 discipline for the NEW (e)/(f) tests

- **Mode-8 (`-O` strips bare `assert`).** EVERY new (e)/(f) tripwire is `-O`-safe:
  `pytest.fail` + `np.testing.assert_array_equal` (function calls, fire under `-O`).
  The §E.1.2 emit-reorder pin, the §E.1.3 kernel-fold sha256 guard, the §E.1.5
  retirement structural assertions, the §E.2.1 grep-gate, and the §E.2.2 audit are
  ALL function-call-based. Grep `^\s*assert ` in every new (e)/(f) test file →
  expect zero. The §E.2.1 collection + `-W` Sphinx gates are not assertions at all
  (collection success / build success) — Mode-8-immune by construction.
- **Mode-9 (FP-invariance only in a degenerate regime).** (e) and (f) are
  RE-LAYERING + RELOCATION, NOT a splitting/accelerator — the converged FP does NOT
  change, and they FORK NO default. The classic Mode-9 splitting trap does not
  directly apply. BUT every (e) equivalence gate (the emit-reorder pin §E.1.2)
  MUST run on an **ANISOTROPIC + HETEROGENEOUS + ≥2G NON-SQUARE** config (12×7 /
  5×9 + synthetic d=3), NEVER the isotropic-reflective square box: the octant
  interior recurrence + the moment emit + the x↔y face mapping ALL degenerate on a
  1G flat square (vv §H1/§H2 + [[phase5a_5c_angular_windowing_lessons]] SQUARE-
  snapshot x↔y blindness + THE TRAP — a P0-isotropic config lets an ℓ≥1-dropping
  emit-collapse pass silently; the FULL moment tensor assertion in §E.1.2 is the
  catcher). The NON-SQUARE shapes are LOAD-BEARING — they catch a mu_x↔mu_y swap
  (vv failure mode #2) the collapsed walk could introduce. The d=3 synthetic shape
  guards against a 2-D hardcode creeping into the parameterized walk (C3 regression
  — the amendment's reason the walk is `_OctantWalk` d-generic, not `_OctantWalk2D`).
