---
name: c3-6-honest-d-dispatch-verification
description: C3.6 honest d-dispatch V&V plan (#220, N-D tail). (a) sweep_schedule d-generic _octant_sweep(entry,ndim)+_outgoing_faces enumerate(signs) closes the d=3 silent z-face shed hole (ERR-056 class) — synthetic-d3 schedule pin via PURE helpers (no SNMesh, entry duck-typed on .label/.indices); (b) SNMesh.streaming(axis) d-generic = bit-id by construction (mu_x IS axis_cosines(0)); (c) ScanMarch.supports narrow is_1d or (is_cartesian and ndim==2) → d=3 falls through to FullFieldWavefront. ⭐ WORKING TREE IS MID-CARVE: (a) prod half-written, callers don't pass ndim (TypeError), 7 tests red from shim removal. Mode-8: schedule tests bare-assert inert under -O.
metadata:
  type: project
---

C3.6 verification plan — **honest d-dispatch** (#220), the N-D campaign tail
after the #222 sweep-strategy carve. Worktree `worktree-sn-nd-layout`, HEAD
`0036acc` (#222 CLOSED). Extends [[c3-dim-generic-sweep-dag-verification]] (the
synthetic-d3 idiom, the x↔y-swap moat), [[sweep-strategy-carve-verification]]
(the `is_cartesian`+`ndim` orthogonality, `default_for` registry order), and
[[scan-march-verification]] (the ScanMarch production-default reality, the
ERR-056 shed-order class). This is the proactive `test-architect` gate before
the (a)/(b)/(c) implementation is finished turn-by-turn.

## ⭐⭐ LOAD-BEARING SURPRISE — the working tree is MID-CARVE (changes the plan)

The (a) PRODUCTION change is HALF-WRITTEN in the uncommitted working tree
(`git status`: `sweep_graph.py` + `sweep_schedule.py` both `M`). Verified
this session by reading the diff:

- ✅ `OctantLabel.sign_x`/`.sign_y`/`.streams_in_2d` shims **already removed**
  (`sweep_graph.py`); a new `AXIS_NAMES = ("x","y","z")` SoT added.
- ✅ `_octant_sweep(entry, ndim)` **already** keeps `signs[:ndim]` zero-padded.
- ✅ `_outgoing_faces` **already** derives `f"{AXIS_NAMES[a]}{'max'|'min'}"`
  from `enumerate(label.signs)`, skipping `s==0`.
- ✅ `_OUT_FACE` hand-listed dict **already** deleted.
- ❌ **The two CALLERS were NOT updated**: `jacobi` (`:99`) and `gauss_seidel`
  (`:119`) still call `_octant_sweep(entry)` with NO `ndim` arg → a
  **`TypeError` at runtime**. The schedule is currently BROKEN.
- ❌ The 4 test consumers in `test_sweep_schedule.py` (`_expected_outgoing`
  helper `:83-89` + `test_gs_box_half_reflective` `:187-189`) + the 2
  consumers in `test_sweep_graph.py` (`:46-47`, `:55-60`) still read
  `.sign_x`/`.sign_y`/`.streams_in_2d` → `AttributeError`.

**Live red count at HEAD (with the dirty tree), `python -O`:**
`tests/sn/sweep/core/ + cartesian_2d/` = **7 failed, 464 passed, 2 skipped,
4 xfailed**. The 7 reds: 4 in `test_sweep_schedule.py` (the `.sign_x` ones),
3 in `test_scan_march_end_to_end.py` + `test_affine_carve_bit_identity.py`
(`si_2d_p1_aniso_het` case + G4.a + G4.b — these go through the broken
schedule on the reflective path).

⚠ **The implementer's FIRST move is to FINISH (a):** pass `ndim` from the two
callers (`_octant_sweep(entry, sn_mesh.ndim)`), then migrate the 6 test
consumers (§3). The gates below assume (a) is being COMPLETED, not started.
This is NOT a "the production change is correct, write gates for it" task —
it is a "the production change is half-applied and currently red, the gates
must (i) drive it green, (ii) add the d=3 coverage that justified the carve."

## ⭐ The carve's WHY — the d=3 silent-wrong hole (DEMONSTRATED this session)

The headline of (a) is NOT a refactor — it is a **correctness surface**. Probe
output (`PYTHONPATH=<worktree>`, pure helpers, no mesh):

```
current _octant_sweep on 3-sign entry (+1,-1,+1) -> label OctantLabel(signs=(1,-1))  [truncates sign_z]
current _outgoing_faces (+1,-1)  -> ('xmax', 'ymin')
=> a +z octant NEVER sheds zmax: at d=3 a reflective z-face gets NO G-S reflect group → WRONG fixed point
```

This is the **ERR-056 class** (`numerical-bug-signatures` Sig 1 / vv Mode 9):
a reflective-BC outflow that is silently dropped converges to a stable-but-
wrong fixed point (global balance still telescopes — vv §H3). The original
ERR-056 was the *shed-ORDER* on a shared face (diagonal cubature); this is the
*shed-EXISTENCE* on the third axis (the z-face never enters the schedule at
all). Both are silent-wrong reflective-BC fixed-point bugs. The (a) change
closes the existence hole; the synthetic-d3 schedule pin (G-a3 below) is the
gate that PROVES it closed.

**Claim-layer declaration (vv §1.5).** Every (a) gate is a **flux-shape /
structural** claim (the schedule is a mesh-time combinatorial object; no
eigenvalue, no MMS). Pillar = **closed-form structural** (the outgoing-face
set and the reflect-group assignment are derivable from the octant signs by
hand — first-principles, structurally independent of the schedule module's
internal derivation, vv L11). No eigenvalue claim rides on (a) directly; the
eigenvalue safety is inherited transitively through the d=2 FP-invariance
end-to-end gates (G4.a/G4.b, already in the suite) staying green.

---

## 1. THE GATE LIST PER CHANGE

### (a) `sweep_schedule.py` d-generic — the ERR-056-class correctness surface

| Gate | Level | File (new/ext) | Marker | Pins |
| ---- | ----- | -------------- | ------ | ---- |
| **G-a0** finish-the-carve | — | (production) | — | callers pass `ndim`; suite no longer TypeErrors. NOT a test — the precondition. |
| **G-a1** d=2 bit-identity (jacobi+G-S structure) | foundation | EXT `test_sweep_schedule.py` (migrate existing) | `foundation` | the 11 existing structure tests stay green AFTER migration; d=2 schedule objects byte-identical (same groups/labels/reflect_faces). |
| **G-a2** ERR-056 shared-face (diagonal cubature) | foundation | EXT `test_sweep_schedule.py` (migrate `test_gs_diagonal_*`) | `foundation` + `catches("ERR-056")` | the lebedev shared-face last-group assignment survives the d-generic rewrite (already exists; migrate `.sign_x` → `.signs`). |
| **G-a3** ⭐ synthetic d=3 schedule pin (NEW) | foundation | NEW `test_sweep_schedule_nd.py` | `foundation` | `_octant_sweep`/`_outgoing_faces` keep all 3 signs at ndim=3; a +z octant sheds `zmax`; a reflective z-face gets a reflect group. THE gate that justifies the carve. |
| **G-a4** d=2 zero-pad (slab quad over 2-D mesh, the FEWER-signs arm) | foundation | NEW `test_sweep_schedule_nd.py` | `foundation` | a quadrature with fewer signs than the mesh has axes zero-pads (sign 0 = no streaming on that axis) — the `a < len(label) else 0` branch. |
| **G-a5** d=2 end-to-end FP-invariance (inherited) | l1 | EXT `test_scan_march_end_to_end.py` (drive green) | `l1` + `catches("ERR-056")` | G4.a/G4.b converged ψ* schedule-invariant — currently RED through the broken schedule; (a)-complete drives them green. NO new test; confirm they flip. |

### (b) `SNMesh.streaming(axis)` genuinely d-generic — bit-identity BY CONSTRUCTION

| Gate | Level | File (new/ext) | Marker | Pins |
| ---- | ----- | -------------- | ------ | ---- |
| **G-b1** stencil value bit-id (d=1) | foundation | EXT `test_sweep_regression.py::test_stencil_values_cartesian` (UNCHANGED) | (class `TestSNMesh`, no marker → foundation-by-default; see §6) | `streaming(0)[n,i] == 2|μ_x[n]|/dx[i]` to rtol 1e-14 — STAYS green; the `axis_cosines(0)`-vs-`mu_x` view identity makes it bit-id. |
| **G-b2** DD-denom equivalence (d=2, lebedev) | foundation | EXT `test_sweep_regression.py::test_stencil_dd_denom_equivalence` (UNCHANGED) | (TestSNMesh) | `Σ_t + streaming(0)[n,i] + streaming(1)[n,j] == old hand denom` rtol 1e-14 — STAYS green. |
| **G-b3** ⭐ per-axis view identity (NEW, the bit-id-by-construction pin) | foundation | NEW small test in `test_sweep_regression.py` or `test_snmesh_consumes_reduced.py` | `foundation` | `np.testing.assert_array_equal(sn.streaming(0), 2|axis_cosines(0)|/dx)` AND `sn.streaming(1) == 2|axis_cosines(1)|/dy` — pins that the new `range(ndim)` tuple build reads `axis_cosines(a)` (a VIEW of `mu_x`), so the value is bit-IDENTICAL to the legacy `(streaming_x, streaming_y)[axis]`, not merely close. |
| **G-b4** non-Cartesian raises (regression) | foundation | EXT/NEW | `foundation` | `streaming(axis)` on a sphere/cyl mesh raises `AttributeError` (the curvature guard, `geometry.py:670`) — pins the Cartesian-only contract survives the d-generic rewrite. |
| **G-b5** axis bounds (regression) | foundation | NEW | `foundation` | `streaming(2)` on a 2-D mesh raises `IndexError` (`0<=axis<ndim`, `:676`); pins the new tuple is indexed in-range. |

### (c) honest `supports` + wording sweep — NO numerical change

| Gate | Level | File (new/ext) | Marker | Pins |
| ---- | ----- | -------------- | ------ | ---- |
| **G-c1** ⭐ ScanMarch d=3 NOT-supported (NEW) | foundation | NEW `test_strategy_d_dispatch.py` OR EXT existing strategy-selection test | `foundation` | `ScanMarch.supports(synthetic_d3_mesh).ok is False` — the narrowed predicate `is_1d or (is_cartesian and ndim==2)` refuses d=3. Needs a d=3 mesh-LIKE object (see §2-c). |
| **G-c2** ScanMarch d=2 + d=1 STILL supported (regression) | foundation | SAME file | `foundation` | `ScanMarch.supports(box_2d).ok` AND `.supports(slab_1d).ok` AND `.supports(cyl_1d).ok` — the narrowing did NOT lose d≤2 / 1-D-any-geom coverage. |
| **G-c3** ⭐ default_for(d=3 Cartesian) → FullFieldWavefront (NEW) | foundation | SAME file | `foundation` | with ScanMarch refusing d=3 and MovingFrontierWindow refusing (`ndim==2`), `default_for` falls through the registry to `FullFieldWavefront` (any-d Cartesian spine, `:988`). The fall-through chain is the (c) headline. |
| **G-c4** default_for(d=2) UNCHANGED → ScanMarch (regression) | foundation | SAME file | `foundation` | d=2 default is STILL ScanMarch (the S6.9 Fork-B2 production default) — (c) did NOT disturb the d=2 selection. |
| **G-c5** bit-identity anchors UNTOUCHED (regression) | l2/l3 | EXT `test_affine_carve_bit_identity.py` + A2D-1 source-hash | `l2`/`l3` (existing) | the affine sha256 goldens stay byte-identical — (c) is `supports`-predicate + docstring/message wording ONLY, NO `default_for` reorder, NO kernel touch. **NO regeneration this time** (task constraint). |
| **G-c6** wording sweep (informational) | — | grep gate | — | the narrowed `supports` reason string + any d=3-referencing docstring read truthfully (no "any-d Cartesian" claim left on ScanMarch). Not a pytest test — a review item. |

---

## 2. THE SYNTHETIC d=3 SCHEDULE PIN DESIGN (deliverable item 2)

**The seam — VERIFIED this session.** `SweepSchedule.jacobi/gauss_seidel(sn_mesh)`
read EXACTLY three things off the mesh (grepped, `sweep_schedule.py`):

1. `sn_mesh.quad.octants` — a tuple of `DiscreteMeasurePartition` (consumed by
   `_octant_sweep`, which reads ONLY `entry.label` (a tuple) + `entry.indices`
   (an int array)).
2. `sn_mesh.trace.layout.faces` — the face-name set (consumed by
   `_reflective_faces`).
3. `getattr(sn_mesh, f"bc_{face}")` per face — the BC strings.
   (plus `sn_mesh.ndim` once the carve is finished, for the `_octant_sweep`
   `ndim` arg.)

`_octant_sweep(entry, ndim)` and `_outgoing_faces(label)` are **PURE
functions** on those duck-typed inputs — NO SNMesh needed.
`DiscreteMeasurePartition` has fields `(label, indices, measure)` but
`_octant_sweep` touches only `.label`/`.indices`, so a `SimpleNamespace` (or
a tiny `@dataclass`) duck-type suffices.

**⭐ RECOMMENDATION (do NOT over-engineer): test the pure helpers DIRECTLY.**
The d=3 schedule structure is fully determined by `_octant_sweep` +
`_outgoing_faces` + the `gauss_seidel` last-group-assignment loop. Two tiers:

**Tier 1 (G-a3/G-a4) — pure-helper unit pins (cheapest, primary).** No mesh,
no measure. Forge `entry = SimpleNamespace(label=(+1,-1,+1), indices=np.array([0,1]))`:

```python
def test_octant_sweep_keeps_all_three_signs_at_ndim3():
    entry = SimpleNamespace(label=(+1, -1, +1), indices=np.array([0, 1, 2]))
    sw = _octant_sweep(entry, ndim=3)
    np.testing.assert_array_equal(sw.label.signs, (+1, -1, +1))   # NOT truncated

def test_outgoing_faces_includes_z_at_ndim3():
    assert _outgoing_faces(OctantLabel((+1, -1, +1))) == ("xmax", "ymin", "zmax")
    # ^ use pytest.fail / np.testing — NOT bare assert (Mode-8; see §6)

def test_octant_sweep_zero_pads_slab_quad_over_2d_mesh():   # G-a4
    entry = SimpleNamespace(label=(+1,), indices=np.array([0]))   # 1-sign quad label
    sw = _octant_sweep(entry, ndim=2)
    np.testing.assert_array_equal(sw.label.signs, (+1, 0))        # zero-padded
```

**Tier 2 (G-a3 full) — a small injectable-schedule pin via a duck-typed mesh.**
The current `jacobi/gauss_seidel` take `sn_mesh` and read the 3+1 attrs above.
The cheapest way to drive the WHOLE `gauss_seidel` last-group-assignment loop
at d=3 WITHOUT an SNMesh is a tiny fake:

```python
@dataclass
class _FakeMesh3D:
    quad: object               # SimpleNamespace(octants=(...,))
    trace: object              # SimpleNamespace(layout=SimpleNamespace(faces=(...)))
    ndim: int = 3
    bc_xmin = bc_xmax = bc_ymin = bc_ymax = "reflective"
    bc_zmin = bc_zmax = "reflective"
```

Then `SweepSchedule.gauss_seidel(_FakeMesh3D(...))` exercises the real
production loop at d=3 and the pin asserts: every reflective face
(`{xmin,xmax,ymin,ymax,zmin,zmax}`) appears in EXACTLY ONE group's
`reflect_faces`, and `zmax`/`zmin` are PRESENT (the hole the carve closes).

⚠ **The ONE refactor I recommend (small, principled).** The fake mesh above
must expose `ndim`, `quad.octants`, `trace.layout.faces`, and `bc_*`. That is
a 4-attr duck-type — acceptable, NO production refactor needed. **Do NOT** add
an "injectable mesh reads" abstraction to production (coding-elegance Pattern 6:
1 consumer, the test — defer). The duck-type lives in the test file. **Prefer
Tier 1** (pure-helper) for G-a3/G-a4 (it pins the actual seam the carve
touched); add Tier 2 ONLY for the full `gauss_seidel` last-group-at-d=3
assignment (the ERR-056-existence pin). Both are `foundation`, both
`np.testing`/`pytest.fail`.

**Structural independence (vv L11) — the hand-derived d=3 oracle.** The pin's
expected face set must be derived FIRST-PRINCIPLES, not from the module's own
`_outgoing_faces`. Re-use the existing `_expected_outgoing` pattern in
`test_sweep_schedule.py` (the independent first-principles map `sign_x>0 →
xmax`), EXTENDED to z: `sign_z>0 → zmax`, `sign_z<0 → zmin`. This independent
map is the oracle `_outgoing_faces` is checked against — exactly the L11
discipline the existing slab tests already use. (The §3 migration makes
`_expected_outgoing` d-generic, so it serves both d=2 and d=3.)

**(c) synthetic d=3 mesh-like for G-c1/G-c3.** `ScanMarch.supports(mesh)` /
`MovingFrontierWindow.supports(mesh)` / `FullFieldWavefront.supports(mesh)`
read ONLY `mesh.is_1d`, `mesh.is_cartesian`, `mesh.ndim` (grepped). A 3-attr
`SimpleNamespace(is_1d=False, is_cartesian=True, ndim=3)` drives all three
`supports` + the `default_for` selection at d=3 with NO mesh. `default_for`
additionally calls `cls(mesh)` on the WINNING strategy — so G-c3 must either
(i) assert at the `supports`-level only, or (ii) stub `cls.__init__`.
**Recommend (i)** — assert the three `supports` predicates AND that
`[c for c in LOSS_REPRESENTATIONS if c.supports(fake3d).ok][0] is
FullFieldWavefront` (this IS `default_for`'s "first that applies" contract,
re-derived structurally-independently from the registry order without invoking
the un-constructible `cls(fake3d)`). This also closes F3 below.

---

## 3. SHIM-RETIREMENT TEST-MIGRATION LIST (deliverable item 3)

The shims (`OctantLabel.sign_x`/`sign_y`/`streams_in_2d`) are ALREADY removed
from production; the migration is the remaining TEST consumers. Exact lines,
what each assertion BECOMES:

| File:line | Current | Becomes |
| --------- | ------- | ------- |
| `test_sweep_schedule.py:83-90` (`_expected_outgoing` helper) | `if label.sign_x > 0: faces.add("xmax") elif label.sign_x < 0: ...` (×2 for x,y) | `for a, s in enumerate(label.signs): if s>0: faces.add(f"{('x','y','z')[a]}max") elif s<0: faces.add(...min)`. Keep it a first-principles MAP (NOT a call to `_outgoing_faces` — L11 independence). Now d-generic → reusable by G-a3's d=3 oracle. |
| `test_sweep_schedule.py:187-189` (`test_gs_box_half_reflective`) | `if label.sign_x > 0: assert g.reflect_faces == ("xmax",) elif label.sign_x < 0: ...` | `if label.signs[0] > 0: ... elif label.signs[0] < 0: ...` (read `.signs[0]` directly). |
| `test_sweep_graph.py:42-47` (`TestOctantLabel::test_valid_signs`) | `assert lab.sign_x == sx; assert lab.sign_y == sy` | `np.testing.assert_array_equal(lab.signs, (sx, sy))` — AND convert to `np.testing` (the class is `@pytest.mark.l0` but uses bare `assert` → Mode-8 inert under `-O`; the migration is the moment to fix it). |
| `test_sweep_graph.py:55-60` (`test_streams_in_2d`) | `assert OctantLabel((+1,+1)).streams_in_2d` ×4 + `assert not ...((0,0)).streams_in_2d` | rename → `test_streams` + `.streams` (the real property; `streams_in_2d` was the deprecated alias). Convert bare `assert` → `if not X: pytest.fail(...)` (fires under `-O`). |
| `test_2d_octant_sweep_equivalence.py:285,701,903` | `sign_x == 0 == sign_y` in COMMENTS/strings only (grepped — no live `.sign_x` attribute access) | docstring/comment text fix only: `signs == (0, 0)` or "all-zero label". NO assertion change. |

**Post-migration grep gate (audit deliverable):** `grep -rn "sign_x\|sign_y\|
streams_in_2d" orpheus/ tests/` returns ZERO outside the `AXIS_NAMES`-derived
`f"{axis}min/max"` strings and the `signs[0]/[1]` index reads. Confirms the
shim is fully retired (coding-elegance aggressive-retirement + the user's
[[retirement-means-test-migration]]).

---

## 4. REGRESSION MATRIX — before/after EACH commit (deliverable item 4)

The carve is naturally three commits (a)/(b)/(c). Run order + expected counts
(`python -O`, sequential — xdist deadlocks per [[sn-taxonomy-reorg-mapping]]):

**Baseline NOW (dirty tree, HEAD `0036acc`):** `tests/sn/sweep/core/ +
cartesian_2d/` = **7 failed, 464 passed, 2 skipped, 4 xfailed** (the 7 reds
are the in-flight shim breakage, §SURPRISE). The smaller schedule+graph+
streaming batch (`test_sweep_schedule.py test_sweep_graph.py
test_snmesh_consumes_reduced.py test_sweep_regression.py`) = **108 passed**
ONLY because the bare-assert schedule tests pass vacuously under `-O`; the 4
`.sign_x`-via-`_expected_outgoing` ones go red the moment `cartesian_2d/` is
co-run (collection of the full path). Treat 108 as the streaming/graph anchor,
the 7-red figure as the schedule reality.

**Commit (a) — finish d-generic schedule + migrate consumers + add d=3 pins:**
- BEFORE: the 7 reds above.
- Order: (1) `test_sweep_schedule.py` (the 4 `.sign_x` reds → green after
  caller fix + consumer migration; +2 NEW d=3 pins G-a3/G-a4 if co-located, or
  in the new `test_sweep_schedule_nd.py`). (2) `test_sweep_graph.py` (2 consumer
  migrations green). (3) `test_scan_march_end_to_end.py` (G4.a/G4.b flip green —
  the schedule un-breaks the reflective path). (4)
  `test_affine_carve_bit_identity.py` (`si_2d_p1_aniso_het` flips green —
  same reflective-schedule path; **byte-identical**, NO regeneration). (5) full
  `tests/sn/sweep/core/ + cartesian_2d/`.
- AFTER (expected): **0 failed**, ~466-470 passed + the NEW d=3 pins (≥3). The
  4 xfailed stay xfailed.
- ⚠ The 4 schedule reds that are bare-`assert`-elsewhere and pass under `-O`
  while asserting nothing: the migration to `np.testing` will make them
  ACTUALLY fire — run `test_sweep_schedule.py` ALSO **without `-O`** once to
  confirm the migrated assertions are live (Mode-8).

**Commit (b) — `SNMesh.streaming(axis)` d-generic:**
- Anchor suites: `tests/sn/sweep/core/test_sweep_regression.py` (G-b1/G-b2,
  the `TestSNMesh` stencil tests — MUST stay green, bit-id by construction),
  `tests/sn/primitives/test_snmesh_consumes_reduced.py` (G-b4, the `streaming_x
  is not None` pins `:174-175`).
- BEFORE/AFTER: **byte-identical** — `mu_x IS axis_cosines(0)` (a `@property`
  view, confirmed `directional.py:294-302`), so the `range(ndim)` tuple build
  produces the SAME arrays. Add G-b3 (the explicit `assert_array_equal` view-
  identity pin) + G-b5 (axis-bounds IndexError).
- Combined `test_sweep_regression.py + test_snmesh_consumes_reduced.py` are part
  of the 108-pass anchor batch (clean lines). After (b): +2 NEW (G-b3/G-b5),
  all green.

**Commit (c) — narrow `ScanMarch.supports` + wording:**
- Anchor suites IN ORDER: (1) NEW `test_strategy_d_dispatch.py` (G-c1..G-c4).
  (2) `test_affine_carve_bit_identity.py` (G-c5 — **byte-identical**, the
  predicate narrowing does NOT change d=2 selection so the production path is
  untouched; **NO golden regeneration** — task constraint). (3) A2D-1 source
  hash `test_streaming_operator.py::TestT4dApply2DCartesianSourceHashPin`
  (untouched — `supports` is not in `_apply_2d_cartesian`'s source). (4)
  `test_scan_march_end_to_end.py` (G6 end-to-end — d=2 default still ScanMarch).
  (5) `tests/sn/eigenvalue/test_keff_2d.py` (the eigenvalue safety leg).
- `test_affine_carve_bit_identity.py + test_scan_march_end_to_end.py +
  test_streaming_operator.py` clean-tree baseline this session = **68 passed**
  + the 3 (a)-red ones; after (a)+(c): all 71 green.

**Final full-campaign ladder (the C3.6 acceptance run):** the fast ladder from
[[c3-dim-generic-sweep-dag-verification]]:
`primitives/ sweep/core/ sweep/slab/ sweep/curvilinear/ solve/` +
`eigenvalue/test_keff_2d.py` + `operators/` — **expect the prior
862-pass/31-xfail baseline + the new C3.6 pins, 0 failed.** Deselect
`test_keff_slab.py::test_heterogeneous_absolute_keff` (#212 hang). Standing
held reds untouched (all `xfail`): #206 cyl-matvec, #195 MMS@160. Run
SEQUENTIAL (xdist 3.8.0 + Py3.14 loadscope deadlock).

---

## 5. FAILURE MODES THE GATES ABOVE DO NOT CATCH — + cheapest closer (item 5)

| # | Foreseen miss | Why the gates above are blind | Cheapest closing gate |
| - | ------------- | ----------------------------- | --------------------- |
| **F1** | ⭐ **Mode-8 inert schedule tests.** `test_sweep_schedule.py` + `test_sweep_graph.py::TestOctantLabel` use BARE `assert` under the canonical `-O` invocation → they assert NOTHING. Confirmed: only 4 of 11 G-S structure tests went red on shim removal (the `.sign_x`-via-`_expected_outgoing` ones raise `AttributeError` at CALL regardless of `-O`); a pure-`assert` structure test that silently passed could be asserting a WRONG d=2 schedule and we'd never know. | The migration (§3) MUST convert the touched bare-`assert` lines to `np.testing.assert_array_equal` / `pytest.fail`. G-discipline: grep `^\s*assert ` in the migrated files; expect zero on the lines the carve touches. Run `test_sweep_schedule.py` once WITHOUT `-O` in the (a) acceptance to confirm the migrated assertions are live. |
| **F2** | **x↔y↔z axis swap in the face-name derivation.** `_outgoing_faces` now indexes `AXIS_NAMES[a]` by `enumerate(label.signs)`. A transposed `AXIS_NAMES = ("y","x","z")` or an off-by-one would produce `ymax` where `xmax` belongs — INVISIBLE on a square/symmetric ALL-reflective box (a swapped name still maps to a reflective face). The synthetic d=3 pin (G-a3) on an all-reflective box is blind to this. | G-a3 must use a **mixed-BC** config: x-reflective, y-vacuum, z-reflective. Then a `(+1,+1,+1)` octant must reflect ONLY `xmax`,`zmax` (y vacuum). An x↔y swap reflects the wrong face → caught. This is the d=3 analog of the existing `test_gs_box_half_reflective` 2-D pin (kept by §3). One mixed-BC synthetic d=3 pin closes it. |
| **F3** | **`default_for` fall-through asserted only at `supports`-level (G-c3 rec (i)).** Asserting the three `supports` predicates does NOT prove `default_for` RETURNS `FullFieldWavefront` at d=3 — a registry-order bug could pass while changing the d=2 default. | Closed BY the G-c3 (i) form itself: assert `[c for c in LOSS_REPRESENTATIONS if c.supports(fake3d).ok][0] is FullFieldWavefront` (re-derives default_for's "first applies" contract on the REAL registry tuple, no construction). Plus G-c4 (d=2 still ScanMarch) guards the other direction. Free — same file. |
| **F4** | **The d=3 schedule pin proves STRUCTURE, never the d=3 SWEEP VALUE.** A correct schedule fed to a wrong d=3 sweep kernel still converges wrong; no d=3 mesh exists, so no d=3 fixed-point gate is possible. The carve is structurally d-generic but the VALUE claim at d=3 is DEFERRED to Mesh3D (C5). | STATE this boundary explicitly in the d=3 pin docstrings (vv Mode 7 — declare what the synthetic config does NOT activate): "G-a3 pins the d=3 schedule COMBINATORICS only; the d=3 converged-value claim is C5/Mesh3D, gated by an end-to-end FP-invariance + k_inf at that phase." No new gate — a documented scope boundary so a future reader does not mistake structure-coverage for value-coverage. (Honest analog of the `ScanMarch` d=3 DEFERRAL — its kernel hard-unpacks `ng,nx,ny`; both (a) and (c) are "d-general structure, d=2 values".) |
| **F5** | **(b) bit-identity claim rests on `mu_x IS axis_cosines(0)` staying a VIEW.** If a future refactor makes `mu_x` a COPY (or `axis_cosines` returns a fresh array with different FP rounding), G-b1/G-b2's rtol=1e-14 would STILL PASS (values math-equal) but the "bit-identical by construction" claim silently degrades to "close". | G-b3 (the NEW `np.testing.assert_array_equal` — EXACT, not `allclose`) is precisely this guard: it FAILS if `streaming(0)` ever stops being bit-identical to the legacy `streaming_x`. Already in the (b) list; flagged so it is not dropped as "redundant with G-b1" — it is NOT (G-b1 is rtol-based, G-b3 is exact). |

---

## 6. VV TAGGING + DISCIPLINE

- **Markers (per [[feedback-vv-tagging]]):** the schedule/streaming/strategy
  structure pins carry NO equation `:label:` → `@pytest.mark.foundation` (the
  schedule is a mesh-time combinatorial object; the streaming accessor is a
  data-structure invariant; `supports`/`default_for` is factory-output). The
  d=2 end-to-end FP-invariance ones STAY `l1` + `catches("ERR-056")`. **NEVER**
  put `verifies(...)` on the foundation pins. (`test_sweep_regression.py`
  `TestSNMesh` carries no explicit marker today — confirm the harness defaults
  it to foundation, or add `@pytest.mark.foundation` on the new G-b3/G-b5.)
- **Mode-8 (`-O`-strip) — the LIVE hazard this carve INHERITS.** The existing
  schedule/graph structure tests are bare-`assert` and inert under the canonical
  `-O`. EVERY new C3.6 pin (G-a3/G-a4/G-b3/G-b5/G-c1..G-c4) MUST use
  `np.testing.assert_*` / `pytest.fail` / `raise AssertionError`. The §3
  migration is the moment to convert the touched legacy bare-asserts too.
- **Self-improvement (fires BEFORE delivery):** the d=3 silent-z-face-shed hole
  is a NEW INSTANCE of the ERR-056 class (silent-wrong reflective-BC fixed
  point) but a DISTINCT mechanism — shed-EXISTENCE (the z-face never enters the
  schedule) vs the original shed-ORDER (shared-face premature reflect). This is
  NOT a new failure-MODE row for `vv-principles` (it is Mode 9 / Sig 1, tabled).
  **No ERR entry**: d=3 is unconstructible today, so NO production bug shipped —
  the hole is **closed-by-construction before d=3 exists**, not caught in
  production (no "log every caught bug" trigger). G-a3 is tagged
  `catches("ERR-056")` because it is the SAME class the existing entry names;
  the existence-variant rides the existing entry. If a future d=3-capable build
  ever ships the drop, THAT becomes a new ERR.

## 7. PRE-READS (file:line, VERIFIED this session)

- `orpheus/sn/sweep_schedule.py:99,119` — the two `_octant_sweep(entry)` callers
  NOT yet passing `ndim` (the broken-tree TypeError site).
- `orpheus/sn/sweep_schedule.py:165` (new `_octant_sweep(entry, ndim)`),
  `:182` (new `_outgoing_faces` `enumerate(label.signs)`),
  `:204-211` (`_reflective_faces` — the `bc_{face}` mesh-read seam).
- `orpheus/sn/sweep_graph.py:88` (new `AXIS_NAMES` SoT), `:100-146`
  (`OctantLabel` with shims REMOVED).
- `orpheus/numerics/quadrature/directional.py:265` (`axis_cosines`), `:294-302`
  (`mu_x` IS `axis_cosines(0)` — the (b) bit-id-by-construction proof).
- `orpheus/numerics/measure.py:728` (`DiscreteMeasurePartition(label, indices,
  measure)` — `_octant_sweep` touches only label/indices → duck-typeable).
- `orpheus/sn/geometry.py:653-680` (`streaming(axis)` — the (b) seam; `:670`
  curvature-guard AttributeError, `:676` axis-bounds IndexError, `:680` the
  `(streaming_x, streaming_y)[axis]` body to replace), `:1337-1347`
  (`_setup_cartesian` hand-build of `streaming_x`/`streaming_y`).
- `orpheus/sn/loss_representation.py:1200-1204` (`ScanMarch.supports` =
  `is_1d or is_cartesian` — the (c) narrowing site), `:1266` (`ng,nx,ny`
  hard-unpack — why d=3 must fall through), `:786-790`
  (`MovingFrontierWindow.supports` = `is_cartesian and ndim==2`), `:984-988`
  (`FullFieldWavefront.supports` = `is_cartesian`, any-d spine),
  `:1436-1441` (`LOSS_REPRESENTATIONS` registry order), `:1444` (`default_for`),
  `:627,2267` (the two `streaming(a) for a in range(ndim)` production consumers).
- Test consumers to migrate: `test_sweep_schedule.py:83-90,187-189`;
  `test_sweep_graph.py:42-47,55-60`; `test_2d_octant_sweep_equivalence.py:
  285,701,903` (comments only). Streaming pins: `test_sweep_regression.py:
  102-136` (`TestSNMesh`), `test_snmesh_consumes_reduced.py:174-175`.
- Existing FP-invariance e2e (drive green): `test_scan_march_end_to_end.py`
  G4.a (`:205`) + G4.b (`:231`, `catches("ERR-056")`).
