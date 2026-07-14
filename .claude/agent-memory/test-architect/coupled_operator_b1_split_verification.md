---
name: coupled-operator-b1-split-verification
description: The B.1a/B.1b gate spec for the Coupled Block Operator campaign (2026-07-10) — the DATA-TYPE-minting precursor to Step-4's operator machinery. B.1a = relax Composite[Interior,Boundary] bounds BulkField/BoundaryField → Field, move the isinstance guards DOWN into FullField.__post_init__ (bit-identical, relocated). B.1b = additive-first split of the unified RadialCharacteristicFlux into RadialCharacteristicInteriorFlux (cells, G=V_cell) ⊕ RadialCharacteristicBoundaryFlux (corner, G=V(R)) + two spaces + two SNMesh props + the named RadialCharacteristicComposite (System B). Three NEW gates: G-A (relaxed-base neutrality, System-A guard relocated message-identical), G-B (THE split-fidelity keystone — arange-fingerprint per-(level,sign) bit-identity + metric-partition + MULTI-LEVEL synthetic crux), G-C (System-B intrinsic composite laws + the 3rd non-Bulk/Boundary multi-instantiation + presence). Sibling of coupled_operator_step4_verification.md (this EXTENDS it for the pre-Step-4 carrier split; N1/N2/M/A2/E4 rows carry over, cited not duplicated).
metadata:
  type: reference
---

# Verifying B.1a / B.1b — the relaxed `Composite` base + the additive ray-block split

**When this applies.** BEFORE Step-4's N-general block machinery, the campaign
first (B.1a) **relaxes** the generic carrier's TypeVar bounds and (B.1b) **splits**
the unified ψ½ ray leaf into a typed interior⊕boundary pair with its own named
composite. This is pure **data-type minting** — NO operator matvec, NO eigenvalue —
so the operator-level rows (Step-4 memo M1-M5, A2a/A2b, E4 two-anchor) stay in the
sibling memo for later sub-steps. This memo designs ONLY the three NEW gates the
Step-4 memo did not anticipate. It EXTENDS `coupled_operator_step4_verification.md`;
its N1 (per-consumer VALUE not proxy) and N2 (multi-instantiation config-blindness)
discipline carry over — cited, not duplicated.

## What B.1a/B.1b actually change (arithmetic-verified against the tree)

- **B.1a** — `orpheus/transport/full_field.py`: `Interior = TypeVar(bound=BulkField)`
  / `Boundary = TypeVar(bound=BoundaryField)` (`:155-158`) relax to `bound=Field`;
  the `isinstance(interior, BulkField)` / `isinstance(boundary, BoundaryField)`
  raises MOVE from `Composite.__post_init__` (`:226-237`) DOWN into
  `FullField.__post_init__` (`:557`, currently `super().__post_init__()` +ψ½). The
  **mesh-identity** check (`:242-251`) STAYS on the base. Message text is unchanged:
  `f"{type(self).__name__}: bulk must be a BulkField"` — `type(self)` already
  resolves to the concrete class, so a FullField still reads `"FullField: bulk must
  be a BulkField"` byte-for-byte. WHY the relaxation is forced: System B's ray leaves
  are `FaceField` siblings (NOT `BulkField`, NOT `BoundaryField`), so the common
  ancestor with a System-A leaf is `Field`.
- **B.1b** — split the unified `RadialCharacteristicFlux(FluxRole,
  RadialCharacteristicField)` into two typed leaves ALONGSIDE it (production
  untouched):
  - `RadialCharacteristicInteriorFlux` — the `cells` legs, `(level,sign)→(ng,nx)`,
    metric `G = tile(V_cell, ng)` per leg (A_BB's marched state).
  - `RadialCharacteristicBoundaryFlux` — the `corner` legs, `(level,sign)→(ng,)` at
    r=R, metric `G = full(ng, V_cell[-1])` per leg (B_b's BC locus).
  - Two spaces (`RadialCharacteristic{Interior,Boundary}Space`) splitting
    `RadialCharacteristicSpace.for_levels` (`radial_characteristic_space.py:251-351`;
    the metric build at `:332-340` — `tile(cell_volumes, ng)` for cells,
    `full(ng, corner_gauge)` for corner, `corner_gauge = cell_volumes[-1]`).
  - Two `SNMesh` cached props mirroring `radial_characteristic_space`
    (`augmented_mesh.py:831-859`).
  - `RadialCharacteristicComposite(Composite[RadialCharacteristicInteriorFlux,
    RadialCharacteristicBoundaryFlux])` — thin subclass: leaf-type `__post_init__`
    guard + a `from_mesh`/`zeros` factory that RAISES on a non-ray-carrying mesh
    (presence = ray-carrying, R12a).

**The unified layout is a DE-INTERLEAVE** (arithmetic-confirmed on the sphere,
ng=2/nx=4/1 level, shape 20): `[−1 cells(8), −1 corner(2), +1 cells(8), +1 corner(2)]`.
The split gathers NON-contiguous slices — `unified[0:8]→interior`, `unified[8:10]→
boundary`, `unified[10:18]→interior`, `unified[18:20]→boundary`. This is exactly
where an offset bug lives, and why G-B's fingerprint is arange-per-slot.

## The gate table (the NOVEL rows — G-A / G-B / G-C)

| Gate | Defining law (claim layer: **structural** — data-type minting) | Oracle (structurally-independent ground) | Tooth (in-process monkeypatch, `-O`) | file:line-estimate | Expected RED |
|---|---|---|---|---|---|
| **G-A1** | System-A guard RELOCATED, message-identical: a non-BulkField interior / non-BoundaryField boundary on a `FullField` still raises the SAME `TypeError` — now from `FullField.__post_init__` | POSITIVE control (valid `FullField(AngularFlux,AngularBoundaryFlux)` constructs) + NEGATIVE VALUE `match="bulk must be a BulkField"` / `"boundary must be a BoundaryField"` (N1: assert the raised MESSAGE, not "some TypeError") | monkeypatch `FullField.__post_init__` → variant dropping the isinstance guards (keeps `super()`+ψ½) → negative stops raising | `full_field.py:557` | negative test constructs instead of raising |
| **G-A2** | Base RELAXATION took effect: `Composite` NOW admits a non-Bulk/non-Boundary `Field` pair | `Composite(interior=RadialCharacteristicFlux, boundary=RadialCharacteristicFlux)` on ONE sphere mesh CONSTRUCTS (both leaves are `Field`, neither Bulk/Boundary) | re-add `isinstance(interior, BulkField)` raise to `Composite.__post_init__` | `full_field.py:226` | the non-Bulk construction REDs `"bulk must be a BulkField"` (proves the check genuinely LEFT the base) |
| **G-A3** | Base RETAINS its Field-level guards (relaxation didn't over-reach) + A2.2 scalar stays green | non-Composite partner `a + 42 → match="same-class partner"`; mesh-mismatch (two distinct meshes) `→ match="share mesh identity"`; **`tests/transport/test_composite.py` stays 0-red** (regression-guard) | monkeypatch `Composite.__post_init__` → drop the base mesh check → the mismatch construction succeeds | `full_field.py:242` | mesh-mismatch stops raising (mesh guard proven still load-bearing on the base) |
| **G-B1** ★ | **THE split-fidelity keystone (forward):** per-(level,sign) `interior.leg == unified.cells` AND `boundary.leg == unified.corner`, bit-identical | unified buffer = `np.arange(space.shape[0])` (L20 fingerprint) → each `(level,sign,part)` slot is a KNOWN unique contiguous range; assert `array_equal` against BOTH the split leaf AND the hand-known arange range | (a) sign-swap; (b) cells↔corner; (c) wrong-offset; (d) **per-level offset `2*pos→pos` / `2*level`** in the split projection/space | split projection + `radial_characteristic_space.py` split-`for_levels` layout | wrong fingerprint (see split-fidelity ledger) |
| **G-B2** | Recompose (reverse) is EXACT: `recompose(interior, boundary).values == unified.values` | `array_equal` against the original unified arange buffer (lossless + invertible) | recompose drops / mis-orders a leg | split recompose fn | round-trip buffer differs O(1) |
| **G-B3** | Metric-partition (ERR-067-adjacent gauge): interior metric per leg == `tile(V_cell, ng)`; boundary == `full(ng, V_cell[-1])`; interleaved sum == unified `inner_product_weights` | HAND-KNOWN formula on the sphere's non-uniform `V_cell = [0.065, 0.458, 1.244, 2.422]` (confirmed distinct) — `tile(V,ng) ≠ full(*, V[-1])`; **assert the weight ARRAY** (Mode-12-safe, not a scalar) | (a) interior metric uses `corner_gauge`; (b) boundary metric uses `V_cell` | `radial_characteristic_space.py:332-340` split | metric weight array differs (V_cell vs V(R)) |
| **G-C1** ★ | **Multi-instantiation crux:** `RadialCharacteristicComposite` is the 3rd `Composite` instantiation, FIRST with non-Bulk/non-Boundary leaves, built from the SAME generic path as System A / A2.2 | `type(sys.interior) is RadialCharacteristicInteriorFlux` AND `type(sys.boundary) is RadialCharacteristicBoundaryFlux` from one generic body; System A `AngularFlux` + A2.2 `ScalarFlux` STAY green | hardcode `AngularFlux` / `isinstance(BulkField)` in the relaxed generic body (`zeros`/`_recombine`/`from_flat` OR `__post_init__`) | `full_field.py` generic body | System B REDs / System A + A2.2 scalar GREEN — the **asymmetry IS the config-blindness evidence** (Step-4 N2, now real domain types) |
| **G-C2** | Presence = ray-carrying (Pattern-4, mesh-keyed, NET-NEW teeth): sphere constructs; cyl-LS / cyl-product / slab RAISE | POSITIVE (sphere → 2 signs × 1 level) + NEGATIVE (each non-carrying builder → `RadialCharacteristicComposite.from_mesh/zeros` raises `match="no starting-direction"`), each with its positive control | drop the presence check from the composite factory | new `RadialCharacteristicComposite.zeros` | a non-carrying mesh constructs (should RAISE) |
| **G-C3** | Intrinsic vector-space laws (mirror `test_composite.py`): affine-reject, scalar scale, flat round-trip, copy, partner reject, mesh | `flux+flux → TypeError match="cannot add two"` (`"affine"` also matches); `2*a`/`a/2`/`-a`; `to_flat`/`from_flat` + wrong-size reject; `copy` value-independence; `a+42 → "same-class partner"`; `c.mesh is mesh` | these ride the BASE algebra — the existing base dunder teeth red here too | base `Composite` dunders | law violated |
| **G-C4** † | Torsor mint (CONDITIONAL — see Displacement flag): `flux−flux → RadialCharacteristicInteriorDisplacement` per block; `flux+Δ → flux`; `a+(b−a)==b` (FP, `atol=1e-15`) | per-block displacement TYPE + values `array_equal`; the torsor LAW to FP round-off | swap the minted sibling / drop the block-thread | split Displacement siblings | wrong displacement type / value |

★ = load-bearing (the split-fidelity keystone G-B1 + the config-blindness crux G-C1).
† = CONDITIONAL on the Displacement split leaves (flag below).

## The split-fidelity mutation ledger (G-B — the additive-split licence)

Config (MANDATORY): **arange fingerprint** on the unified buffer (each slot a unique
range) + **two fixtures** — the sphere-GL production mesh (1 level) AND a **MULTI-LEVEL
synthetic space** `for_levels((0, 2, 5), ng=2, nx=4, cell_volumes=graded)` bypassing the
mesh (non-contiguous level values distinguish enumeration-index `pos` from level-value).
`nx=4 ⟹ cells(ng,4)=8 ≠ corner(ng)=2`, so cells↔corner is a loud size mismatch too.

| tooth | what it mutates | single-level sphere | multi-level synthetic | why the config catches it |
|---|---|---|---|---|
| **sign-swap** | split reads `sign=+1` into `sign=−1`'s slot | RED | RED | arange gives distinct per-sign data (flat data would be BLIND — §0.6) |
| **cells↔corner** | split maps `part="corner"` to the interior leaf | RED (size 8≠2) | RED | nx=4 makes the confusion a shape mismatch (loud) |
| **wrong-offset** | de-interleave reads the wrong contiguous slice | RED | RED | arange fingerprint pins each slot's exact range |
| **per-level offset `2*pos→pos`** | split space / projection level-offset term | **GREEN (blind!)** | **RED** | the sphere carries 1 level → `pos∈{0}` → `2·0=0` DEAD; the asymmetry IS the L20/L11 config-blindness evidence |
| **`2*level` vs `2*pos`** | offset keyed on level-VALUE not enumeration-index | GREEN | RED (level 5 lands at 10·per_sign) | non-contiguous `(0,2,5)` distinguishes value from index |
| **metric V_cell↔V(R)** (G-B3) | interior gets `corner_gauge` / boundary gets `V_cell` | RED (sphere V_cell non-uniform) | RED | sphere shells make `tile(V,ng) ≠ full(*,V[-1])` bit-exactly |
| **recompose drop-leg** (G-B2) | reverse omits/mis-orders a leg | RED | RED | round-trip vs the original arange buffer |

**Why G-B is NOT a self-referential tautology.** The comparison is not merely
"split == unified" (which would inherit the unified's bugs). The oracle is TWO
structurally-independent hand-known grounds: (1) the **arange fingerprint** — I know
`unified.cells(0,-1) == arange[0:8]` by construction, so the split is checked against a
value I computed by hand, not against the unified's own reader; (2) the **`tile`/`full`
formula** on the sphere's independently-known non-uniform `V_cell`. The unified field's
OWN correctness (V_cell metric, layout) is already anchored by the existing carrier
gates (`test_radial_characteristic_carrier.py::TestA4…`) + the ERR-067 Mode-12 closure,
so the chain [unified anchored] → [split ≡ arange/formula (G-B)] → [split verified] has
no gap. This is bit-identity by INHERITANCE made sufficient by an independent fingerprint
(the Step-4 E4 discipline, adapted: no operator anchor needed for pure data-type minting).

## Refutation — the additive-split-is-not-a-twin argument (is there a silent divergence G-B misses?)

1. **THE MISS if you only test the sphere: the per-level offset term is DEAD at 1
   level.** The sphere carries exactly ONE seed level, so any `2*pos`/`2*level`
   level-offset term in the split space's layout or the de-interleave projection is
   never activated (`pos∈{0}`). A System-B-carrying-mesh-only suite is BLIND to it —
   the exact L20 `2*pos→pos` config-blindness. **Resolution:** the MULTI-LEVEL synthetic
   fixture `for_levels((0,2,5),…)` (mesh-bypassed, as L20 did) is MANDATORY, and the
   `2*pos→pos` mutation must RED it while staying GREEN on the sphere — the asymmetry is
   the evidence, not a nuisance. Without it, "split IS the unified re-typed" is a
   1-level proxy (Mode-12/L11). **This is the single highest-value refutation.**
2. **Flat/uniform data nulls the sign-swap + de-interleave teeth.** If the unified
   buffer were uniform (or the sphere's `V_cell` accidentally uniform), a sign-swap or
   an offset shuffle would be INVISIBLE (§0.6). Resolution: arange fingerprint (distinct
   per slot) for the layout; the sphere's naturally non-uniform shell volumes for the
   metric (confirmed distinct — no manufactured graded mesh needed on the sphere; but
   the synthetic multi-level fixture MUST pass a graded `cell_volumes` so the metric
   partition is non-degenerate there too).
3. **The split SPACES re-derive the layout+metric — a twin-path hazard.** The two new
   `for_levels` are NEW code that re-computes offsets + the `tile`/`full` metric. If
   they hand-roll offsets (not routed through `FaceLayout.from_named_shapes`) OR re-order
   the named_shapes, they can silently diverge from the unified single-source. G-B1
   (arange, multi-level) catches a layout/order drift; G-B3's **conservation leg**
   (interleaved interior⊕boundary metric == unified `inner_product_weights`) ties the
   split metric back to the unified single-source so a divergent re-derivation reds.
   *Elegance note (Pattern 2):* if the implementer routes the split `for_levels` through
   the SAME `FaceLayout.from_named_shapes` walk (one shared source of the leg order), the
   twin-path is unspellable — flag this as the preferred construction; G-B is the safety
   net if they don't.
4. **cells↔corner silent-equal at nx=1.** At nx=1, `cells(ng,1)` and `corner(ng,)` have
   equal size (ng) → a confusion could be silently equal. Resolution: `nx≥2` (the sphere
   fixture is nx=4) makes the shapes differ — a confusion is a loud mismatch.
5. **What G-B legitimately does NOT cover (out of scope, flag forward):** future
   PRODUCTION divergence — after B.1b the production path still uses the unified field;
   the split is only reached via the projection. If a LATER step populates the split
   leaves INDEPENDENTLY (not via the projection), split and unified could drift. That is
   the Phase-B migration's gate (Step-4 memo P3 "re-type does not move the VALUE" +
   the assemble≡probe keystone), NOT B.1b's. B.1b's contract is exactly "the split IS
   the unified, re-typed" — the retirement licence, nothing more.

## Displacement-sibling flag (G-C4 prerequisite — a concrete design ruling)

`flux − flux → Δ` mints via `Displacement.sibling_of(type(flux))`, keyed by the leaf's
Field-family Rep (`_carrier_rep` = the ONE `Field`-family base, role mixin excluded;
`_displacement.py:171`). Verified mechanics:

- The unified `RadialCharacteristicFlux` has Rep `RadialCharacteristicField`, paired with
  `RadialCharacteristicDisplacement` (same Rep) — that is why the carrier torsor test
  (`test_radial_characteristic_carrier.py::…flux_torsor_triple…`) mints cleanly today.
- The split flux leaves **MUST carry DISTINCT Reps** (new `RadialCharacteristic{Interior,
  Boundary}Field` `FaceField` ABCs, or the concrete leaf is its own Rep). Reusing the
  unified `RadialCharacteristicField` Rep for BOTH split leaves would (a) collide them
  onto ONE `_BY_REP` key and (b) return the UNIFIED-space displacement (wrong space) —
  a real bug the design must avoid.
- For `flux−flux` to mint at all, `RadialCharacteristic{Interior,Boundary}Displacement`
  siblings must be minted AND eagerly imported (the `displacements/__init__.py` populates
  `_BY_REP`). **Absent them, `sibling_of` raises `KeyError`** — the `flux−flux→Δ` and
  `flux+Δ→flux` torsor legs are UNAVAILABLE.

**Ruling (recommend to the main agent):** mint the two Displacement siblings IN B.1b.
The torsor algebra FORCES them (coding-elegance anti-pattern #18; the unified
`RadialCharacteristicDisplacement`'s own docstring says its "existence is FORCED by the
composite torsor algebra" — the same force applies per-block). "No consumer yet ≠
speculative" ([[feedback-defer-only-when-architecture-vague]]): System B iterates the
moment it wires, and the increment IS a per-block displacement. Then G-C4 gates the full
torsor triple. **If deferred:** G-C runs the CORE set (G-C1/2/3 — affine-reject, scalar,
flat, copy, partner, presence, multi-instantiation, ALL displacement-free) and G-C4
becomes a B.1b-followup gate with the gap documented (the "interior-only" fallback the
brief anticipated = the displacement-free CORE, which already fully exercises the
relaxed-base multi-instantiation crux).

## Homes · config · acceptance

- **G-A1** → extend `tests/transport/test_full_field.py` (System-A guard). **G-A2/G-A3**
  → extend `tests/transport/test_composite.py` (base relaxation + retained guards).
  `pytestmark = pytest.mark.foundation` (both files already are).
- **G-B + G-C** → NEW `tests/transport/test_radial_characteristic_composite.py`
  (`@pytest.mark.foundation`). Sphere-GL S4 fixture imported from
  `tests/sn/_test_helpers.py` (`_sphere`, `_cyl_level_symmetric`, `_cyl_product`,
  `_slab`) exactly as `test_radial_characteristic_carrier.py` does; the multi-level
  synthetic space built directly via the split `for_levels`. (Alternative co-location
  `tests/sn/operators/` if the main agent prefers SN-side; the carrier lives in
  `orpheus/transport/`, so transport mirrors the code home — brief's proposal kept.)
- **Config invariants:** carrying = sphere-GL S4 ONLY; cyl/slab = non-carrying CONTROL
  (must RAISE). **≥2G on every value row** (ng=2). Canonical `-O`;
  `np.testing.assert_*` / `pytest.fail` / `pytest.raises`, **never bare assert**
  (Mode-8). Every tooth mutation-verified in-process via `monkeypatch` (NEVER
  `git checkout` an uncommitted file — process-discipline).
- **Blast-radius audits (do before landing):** (a) B.1a — grep every direct
  `Composite(` construction that relied on the BASE rejecting a non-Bulk leaf (only
  tests should; the base is A2.2-new — low risk, but confirm). (b) B.1b — the MRO/
  isinstance blast radius is INHERITED from Step-4 N4 (the `isinstance(_, FullField)`
  consumers): System B is NOT a `FullField`, so leaf-type guards discriminate — the
  new `RadialCharacteristicComposite` must FAIL `isinstance(_, FullField)` and its ray
  leaves must FAIL `isinstance(_, BoundaryField)` (the existing `FaceField`-sibling
  discipline, `_bases.py:1115`).
- **End-to-end acceptance:** `tests/transport` + `tests/sn -m "not slow"` 0-red +
  ratchet transport:1 + `sphinx -W`. The load-bearing deliverables are **G-B1
  (split-fidelity, arange + multi-level crux)** and **G-C1 (the 3rd non-Bulk/Boundary
  multi-instantiation)** — together they LICENSE the additive split + the unified
  field's eventual retirement.
