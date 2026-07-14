---
name: coupled-operator-b2c-builder-verification
description: The DELTA gate spec for campaign sub-step B.2c (2026-07-10) — c1 re-types A_BB (four action surfaces) + A_AB (domain→composite, codomain→full_field_space, emits FullField) onto the composite; c2 mints build_coupled_system → (CoupledOperator, CoupledSpace) 2×2 grid, P1 alignment + P2 presence + the grid≡fused centrepiece + M2-on-real coextensiveness + a real-member block-.H reciprocity. TWO findings against the settled design: (F1) _full_loss_case OMITS B_b_adapter (test_g_adjoint_reciprocity.py:528-542) — the centrepiece must build its OWN complete fused loss WITH B_b on a REFLECTIVE sphere; (F2) radial_characteristic_space == radial_characteristic_composite_space under FunctionSpace.__eq__ (both name="radial_characteristic", shape=(ni+nb,)) → the grid's per-block space check is Mode-12-BLIND to the unified↔composite leaf swap, so "grid constructs" is necessary-NOT-sufficient; the RUNTIME grid.apply + centrepiece is the sufficient catcher. RULING: c1 = array_equal (pure re-label, B.2b §0), but c2 centrepiece bulk rows = principled-equiv rtol=1e-11 (the block-decomposition M-M split, 5.5e-16 — B.2c is NOT uniformly array_equal, unlike B.2b). Extends coupled_operator_step4_verification (4d.2) + coupled_operator_b2b_retype_verification (§0). Sibling of b1/b2b.
metadata:
  type: reference
---

# Verifying B.2c — the co-producing builder + the c1 grid-entry re-types (DELTA)

**When this recipe applies.** A coupled-block campaign has (B.1) posed System B as an
independent `RadialCharacteristicComposite`, (B.2a) landed the semantics-agnostic
`CoupledField`/`CoupledSpace`/`CoupledOperator` machinery, and (B.2b) re-typed `A_BA`
(codomain) + `B_b` (domain/codomain) onto that composite behind transient FullField
adapters. **B.2c now (c1) re-types the two REMAINING System-B-touching blocks — `A_BB` (all
four action surfaces) and `A_AB` (domain + codomain) — and (c2) mints
`build_coupled_system(sn_mesh, mat_xs) → (CoupledOperator, CoupledSpace)`, the typed 2×2
grid.** This is the step-4 memo's **4d.2** (the co-producing builder; SUBSUMES step 6 / task
#34 presence collapse for the grid arm). Governing memos: `coupled_operator_step4_verification.md`
(4d.2 + P1/P2 + the Mode-12 reciprocity deep-dive) and `coupled_operator_b2b_retype_verification.md`
(§0 pure-re-label ruling). **NOTE the step-4 4d.2 line "the ray LEAVES System A's FullField"
is SUPERSEDED** by the B.2-sequencing ruling — the ray eviction is **B.2d**, not B.2c;
through B.2c System A's `full_field_space` stays **3-block** (interior ⊕ trace ⊕
radial_characteristic) on a carrying mesh.

**Verified against HEAD `81b149a`.** All file:line below are read off that tree.

## The object (settled design, arithmetic-verified against the code)

**c1 — grid-entry re-types** (`orpheus/sn/operators/radial_characteristic.py`; engines stay
UNIFIED behind `to_unified`/`from_unified`, the B.2 demote ruling — full unified-leaf
retirement is Phase C/4e):
- **`RadialCharacteristicOperator` (A_BB, :212)** — currently
  `LinearOperator["RadialCharacteristicField"]`, domain=codomain=`radial_characteristic_space`
  (UNIFIED; :305/:309). RE-TYPE **all four action surfaces** — `apply` (:331), `apply_transpose`
  (:368), `solve` (:414), `solve_transpose` (:481) — onto `RadialCharacteristicComposite`
  (domain=codomain=`radial_characteristic_composite_space`), bridging internally
  (`to_unified` → the SAME single-sourced engine → `from_unified`). `inverse()` (:392) body is
  UNCHANGED (`InverseOperator(self)` delegates `.apply`→`solve`/`.solve`→`apply`, so it speaks
  composite automatically). `is_invertible`/`is_adjointable` stay `True`.
- **`RadialCharacteristicSeeding` (A_AB, :580)** — currently
  `LinearOperator["RadialCharacteristicField", "AngularField"]`, domain=`radial_characteristic_space`
  (:723), codomain=`AngularSourceSink._space_for_mesh` (the BULK space, :727). RE-TYPE: domain
  → composite space (parse + `to_unified`, engine unchanged); codomain → **`full_field_space`**
  — `apply` (:734) emits `FullField(interior=<the AngularSourceSink>,
  boundary=AngularBoundarySourceSink.zeros_on, radial_characteristic=<present-zero
  RadialCharacteristicSourceSink>)` (the transitional 3-block presence convention; flips to
  ray=None at B.2d); `apply_transpose` (:792) takes the FullField cotangent, reads `.interior`
  (trace/ray zero-coupled, discarded), returns the composite via `from_unified`.

**c2 — the co-producing builder** (NEW `orpheus/sn/coupled_system.py`, mirroring
`orpheus/numerics/coupled_system.py`):
`build_coupled_system(sn_mesh, mat_xs, *, scattering_order=0) → (CoupledOperator, CoupledSpace)`.
Blocks (loss-sign IN the grid — the block matvec IS the within-group loss action, signs
arithmetic-verified below):
- **(A,A)** = `LC − S − B_a`: `LC = StreamingOperator(sn_mesh) + MultiplicationOperator(coefficient=mat_xs.total_cross_section_field, space=full_field_space)`;
  `S = ScatteringOperator.from_solver_data(mat_xs=, quadrature=sn_mesh.quad, scattering_order=, full_field_space=full_field_space)` (solver.py:917-922);
  `B_a = SNBoundaryOperator(sn_mesh)` (domain=codomain=full_field_space, boundary.py:233/243). An OperatorSum; domain=codomain=full_field_space.
- **(A,B)** = `RadialCharacteristicSeeding(sn_mesh)` **POSITIVE** (its `apply` returns
  `−upstream_numer/V` internally, :787 — it carries the loss-sign). At (0,1): domain=composite,
  codomain=full_field_space.
- **(B,A)** = `−RadialCharacteristicEmission(sn_mesh, S.isotropic_kernel)` — NEGATED via
  `__neg__` → `ScaledOperator(−1.0, Emission)` (operator.py:714; ScaledOperator preserves
  domain=Emission.domain=full_field_space, codomain=Emission.codomain=composite, :1773/:1777;
  is_adjointable inherited True). At (1,0).
- **(B,B)** = `RadialCharacteristicOperator(sn_mesh, mat_xs.total_cross_section_field) −
  RadialCharacteristicBoundaryOperator(sn_mesh)` (`A_BB − B_b`; after c1, A_BB domain=codomain=composite;
  B_b domain=codomain=composite, boundary.py:579/587). An OperatorSum; domain=codomain=composite.
- **Spaces:** `CoupledSpace.from_systems((full_field_space, radial_characteristic_composite_space))`; square (domain=codomain).
- **P2:** `sn_mesh.radial_characteristic_space is None` → 1×1 `[[A_AA]]` over
  `CoupledSpace.from_systems((full_field_space,))`. The block ctors already raise on seedless
  meshes, so System B is simply never built there.

**Config (MANDATORY):** carrying = **sphere-GL S4 ONLY**; cyl/slab = non-carrying CONTROL
(`_cyl_product`/`_slab` imported from `test_loss_transpose_solve`). **≥2G every value row**
(Cardinal Rule). Canonical `.venv/bin/python -O -m pytest -p no:xdist --timeout=300
-p no:cacheprovider` SERIAL; every gate fires under `-O` (`np.testing.assert_*`/`pytest.fail`,
never bare `assert`); every mutation is in-process `monkeypatch` (NEVER `git checkout` an
uncommitted file — L4).

---

## §0 — The claim layer + the ONE ruling that governs the rows (CONTRAST B.2b)

**B.2c is NOT uniformly `array_equal` — this is the sharp difference from B.2b.**

- **c1 re-type value rows = PURE RE-LABEL → `array_equal`** (bit-id INHERITANCE from the
  untouched engine; the bridge `to_unified`/`from_unified` is bitwise, B.1d licence). Identical
  to the B.2b §0 ruling — a `rtol`/`nulp` in a c1 value gate is a RED FLAG.
- **c2's grid≡fused CENTREPIECE bulk rows = PRINCIPLED-EQUIV `rtol=1e-11`** (the campaign's
  measured 5.5e-16 extract-equivalence). The grid `apply` **SPLITS** the fused Morel-Montry
  angular state (`precompute_psi_state` with bulk+seed live) into TWO isolated evaluations —
  seed-zeroed `A_AA` (LC on ψ_A with a present-zero ray) + bulk-zeroed `A_AB` (on ψ_B) —
  summed AFTER. The M-M recurrence is linear in (bulk, seed) jointly, so `f(b,s) = f(b,0) +
  f(0,s)` in EXACT arithmetic; FP reassociation ⟹ principled-equiv. **This is the FIRST
  principled-equiv row in the ψ½ campaign's block realization, and it is INTRINSIC to the
  decomposition** — do NOT try to make it `array_equal` (a tight-`array_equal` centrepiece
  would falsely red the honest block split). The grid IS the "extract" the campaign's "extract
  is principled-equivalent, not bit-identical (5.5e-16)" fact (plan :56-57) describes.

**Claim-layer / pillar (gate on `vv-principles`).** c1 = data-type re-label (structural + bit-id
value). c2-centrepiece = a **construction-equivalence** (grid vs the production fused loss);
its FORWARD-value grounding is INHERITED (the fused loss IS the production loss, already
end-to-end anchored — φ=Q/Σ_t, k_inf — by the existing `tests/sn` wall). The grid cannot SOLVE
at B.2c (`is_invertible=False`; the block solve is step-5), so **no new eigenvalue/fixed-source
anchor is proven here** — E4 solve-anchors DEFER to B.2d (RULING R1). The centrepiece's
forward-equivalence to the anchored fused loss is the B.2c structural tie.

---

## FINDINGS AGAINST THE SETTLED DESIGN (fire before the ink dries)

### F1 — `_full_loss_case` OMITS `B_b_adapter`; the centrepiece must build its OWN complete fused loss

The brief claims the fused reference is "`LC − S − a_ba_adapter − (B_a + B_b_adapter)` — exactly
`_full_loss_case`'s spelling." **WRONG.** `tests/sn/operators/test_g_adjoint_reciprocity.py::_full_loss_case`
(:528-542) spells **`(L+C) − S − SNBoundaryOperator(sn) − a_ba_adapter`** — `B_a` ALONE, no
`B_b_adapter`. (It doesn't need B_b: reciprocity `⟨Aψ,φ⟩_G = ⟨ψ,A.Hφ⟩_G` holds for WHATEVER A
is posed — its own NOTE R2 at :495.) But the grid's **(B,B) = `A_BB − B_b`**, so the fused ray
rows must carry `− B_b` or the centrepiece's ray rows will differ by the whole B_b term.

**RULING:** the centrepiece builds its OWN complete fused loss `fused = LC − S − a_ba_adapter −
B_a − B_b_adapter` (a local helper, NOT `_full_loss_case`), and runs on a **REFLECTIVE** sphere
(`_sphere(bc="reflective")`) so B_b is NON-NULL — on vacuum `_reflect_corner` returns zeros
(`_sphere` docstring :135), which would silently mask the omission. `_full_loss_case` STAYS
untouched (it is the reciprocity fixture, not the centrepiece reference).

### F2 — `radial_characteristic_space == radial_characteristic_composite_space`: the block space-check is Mode-12-BLIND to the unified↔composite leaf swap

The unified `radial_characteristic_space` (name `"radial_characteristic"`,
`radial_characteristic_space.py:417`) and the composite `radial_characteristic_composite_space`
(name `"radial_characteristic"`, `augmented_mesh.py:968`) have **IDENTICAL** `(name,
shape=(n_interior + n_corner,))` — the unified de-interleave and the composite carry the SAME
DOF count. `FunctionSpace.__eq__` is `name == name and shape == shape` (`space.py:170`). So
**the two spaces are `==`**, and `CoupledOperator.__init__`'s per-block check
(`block.domain != domain_members[j]`, `coupled_system.py:650/660`) **PASSES even if A_BB stays
typed on the UNIFIED space** — construction is a **FALSE GREEN**. The runtime, however, CRASHES:
A_BB-on-unified reads `x.values`/`x.cells_view` (the unified layout) but the grid passes it a
`RadialCharacteristicComposite` System-B member.

**Consequences (headline for gate design):**
1. **A_BB's re-type is REQUIRED for RUNTIME correctness, not for construction type-checking**
   (which the `==` masks). The load-bearing catcher is the RUNTIME `grid.apply` on a real
   `CoupledField` (G-c2.1) + the centrepiece (G-c2.3) — NOT "grid constructs" (G-c2.1's
   construction half is necessary-NOT-sufficient, Mode-12-blind).
2. The c1 A_BB re-type value gate (G-c1.1, `array_equal` on the bridged surfaces fed a real
   composite input) IS the proof A_BB speaks composite at runtime.
3. The M2 coextensiveness gate (G-c2.4) is ALSO blind to a unified-for-composite swap (same
   flat size) — it catches a dropped block / wrong layout, NOT the leaf-type. Do NOT credit it
   as the leaf-type catcher.

### F3 — the brief's `is_assemblable` grep-claim is imprecise (conclusion HOLDS)

The brief: "NONE of the real blocks override `is_assemblable`." **Imprecise** —
`MultiplicationOperator` (`multiplication_operator.py:290`, = **C**) and
`IsotropicScattering`/`IsotropicN2N` (`isotropic_scattering.py:301/410`) DO override it. But the
**operative conclusion is CORRECT**: `grid.assemble()` raises `MissingAssembly` because
`StreamingOperator` (L), `ScatteringOperator` (S), and the three ψ½ blocks (A_AB/A_BA/A_BB) do
NOT assemble (base `is_assemblable=False`, `operator.py:655`), so `A_AA = L+C−S−B_a` isn't
(`OperatorSum.is_assemblable = a.is_assemblable and b.is_assemblable`), so the grid isn't, so
`assemble()` raises at block (0,0). Design the M2-on-real re-pin around `to_flat`↔`system_slices`
(as the brief says) — correct target, imprecise premise. ADD a gate asserting
`grid.is_assemblable is False` + `grid.assemble()` raises (documents the unavailability, R2).

### F4 — the A_BB/A_AB re-point blast radius is BROADER than "test_ray_operator — 1 site"

The brief scopes the A_BB re-point to `test_ray_operator.py` (1 construction site, :123). But
**`test_psi_half_coupling.py::TestA_BB_RadialBVP` has ~15 A_BB construction sites** (:807–:1164)
and **`TestA_AB_SeedInjection` ~10 A_AB sites** (:2350–:2579), ALL re-pointing to composite
I/O. The step-4 memo's "TestA_AB_SeedInjection unaffected" note was scoped to the 4d.0 CARRIER
generalization — it does NOT cover the B.2c A_AB RE-TYPE. (F5, minor: the campaign NEXT note's
"A_BB grid entry is direct / A_AB's re-type is a leaf swap" undersells both — A_BB needs the
composite re-type per F2, and A_AB needs a domain swap PLUS a codomain re-type bulk→full_field_space.)

---

## c1 gates — the grid-entry re-types (bit-id INHERITANCE, `array_equal`)

Extend `tests/sn/operators/test_psi_half_coupling.py::TestA_BB_RadialBVP` and
`::TestA_AB_SeedInjection` (the re-points below ARE these gates, re-expressed on composite I/O).

- **G-c1.1 — A_BB four-surface re-type ≡ old unified value.** For each of
  `apply`/`apply_transpose`/`solve`/`solve_transpose`, feed a REAL `RadialCharacteristicComposite`
  (via `_ray_composite`), and assert (i) VALUE `out.to_unified().values == <old-unified-surface>(unified_in).values`
  (`array_equal`); (ii) CONTAINER `type(out) is RadialCharacteristicComposite` (source composite
  for apply/apply_transpose/solve_transpose; flux composite for solve); (iii) `op.domain ==
  op.codomain == radial_characteristic_composite_space`. **This IS the F2 runtime catcher** — a
  still-unified-typed A_BB crashes on the composite input.
  - Teeth: (a) drop the `to_unified` bridge in `apply`/`solve` (return the unified leaf) →
    container assert reds; (b) swap the composite interior↔boundary in `from_unified` → value reds.
- **G-c1.2 — A_BB involution web through composite.** `inverse().apply(composite_src).to_unified()
  == solve(...).to_unified()` and `inverse().solve(composite_flux).to_unified() ==
  apply(...).to_unified()` (`array_equal`); `is_invertible`/`is_adjointable` stay `True`.
- **G-c1.3 — A_AB domain+codomain re-type ≡ old value.** `apply`: `A_AB.apply(composite_seed)`
  → `type is FullField`, `.interior.values == <old>A_AB.apply(unified_seed).values`
  (`array_equal`), `.boundary` present-zero, `.radial_characteristic` present-zero;
  `codomain == full_field_space`, `domain == composite`. `apply_transpose`:
  `A_AB.apply_transpose(FullField_cot).to_unified().values == <old>A_AB.apply_transpose(bulk_cot).values`
  (`array_equal`).
  - Teeth: emit the bulk `AngularSourceSink` directly (drop the FullField wrapper) → `codomain`
    != full_field_space / container reds.

## c2 gates — the builder

Home: NEW `tests/sn/operators/test_psi_half_coupling.py::TestCoupledBuilder` (the step-4 4d.2 home).

- **G-c2.1 — P1 alignment by construction + the RUNTIME proof (F2).** `op, space =
  build_coupled_system(sn, mat_xs)` on the carrying sphere: `space` is a `CoupledSpace` over
  `(full_field_space, composite_space)`; `op.domain == op.codomain == space`; each block at
  `(i,j)` declares the matching member space (assert `op.blocks[0][0].domain == full_field_space`
  DIRECTLY — do NOT rely on the grid's skip-on-None); and **`op.apply(space-shaped
  CoupledField)` RUNS without error** (the sufficient runtime proof; construction alone is
  Mode-12-blind, F2). The builder is the ONLY constructor (mismatched op/space unconstructable).
- **G-c2.2 — P2 presence-structural (NET-NEW teeth).** carrying → `op.n_rows == op.n_cols == 2`
  (System B exists); non-carrying `_cyl_product`/`_slab` → 1×1 `[[A_AA]]`. NEGATIVE:
  monkeypatch the builder's presence predicate to always-build-System-B → on a slab the block
  ctor raises `match=` its specific message (pull EXACT strings): A_BB/A_AB `"carries no
  starting-direction ray (radial_characteristic_space is None)"` (radial_characteristic.py:254/689);
  A_BA `"carries no radial-characteristic ray"` (:1203); B_b `"carries no ψ½ ray
  (radial_characteristic_composite_space is None)"` (boundary.py:559). Each with its positive
  control. (A bare `pytest.raises` false-greens on a downstream crash — `match=` the ctor's own
  message.)
- **G-c2.3 — THE centrepiece: grid ≡ complete fused loss (F1).** Build `fused = LC − S −
  a_ba_adapter − B_a − B_b_adapter` (a local helper WITH B_b, `_sphere(bc="reflective")`, ≥2G).
  Random ψ_full (ray LIVE, via `_random_composite`). Map `coupled = CoupledField(systems=(ψ_A =
  ψ_full with radial_characteristic present-ZEROED, ψ_B =
  RadialCharacteristicComposite.from_unified(ψ_full.radial_characteristic)))`. Assert
  `grid.apply(coupled)` ≡ `fused.apply(ψ_full)` mapped back, **PER-ROW tolerances:**

  | row | claim | bar | why |
  |---|---|---|---|
  | `y_A.interior` | principled-equiv | `assert_allclose(rtol=1e-11)` | the M-M split (seed-zeroed A_AA + bulk-zeroed A_AB summed after); 5.5e-16 measured |
  | `y_A.boundary` | add-order only | `assert_array_almost_equal_nulp(nulp≤8)` / rtol=1e-12 | seed doesn't reach the trace; OperatorSum tree differs by present-zero inserts |
  | `y_A.radial_characteristic` | EXACT zero | `array_equal(zeros)` | A_BB(0)=0, S/B_a/A_AB emit present-zero ray — 0−0−0+0 |
  | `y_B.to_unified()` | add-order only | `assert_array_almost_equal_nulp(nulp≤8)` / rtol=1e-12 | SAME single-sourced bodies on identical inputs (bridge bitwise); `−Emission + (A_BB − B_b)` vs `A_BB − Emission − B_b` |

  Signs arithmetic-verified: bulk `[A_AA_bulk + A_AB] − S − B_a`; ray `A_BB − Emission − B_b`.
  - Teeth: (a) swap A_AB↔A_BA grid placement (offset) → O(1) centrepiece diff (asymmetric); (b)
    drop `− B_b` from (B,B) → ray rows red on the reflective sphere (would GREEN on vacuum —
    F1); (c) the interior-`array_equal` trap: a tightened `array_equal` on `y_A.interior`
    falsely reds the honest M-M split (do NOT ship it).
- **G-c2.4 — M2-on-real coextensiveness + the assemble-arm unavailability (R2).** On the REAL
  grid: per system i, `coupled_field.systems[i].to_flat().size == int(prod(space.systems[i].shape))`
  (FullField 3-block ⊕ composite 2-block); `space.system_slices` extracts + `CoupledField.from_flat`
  round-trips (`array_equal`). PLUS `grid.is_assemblable is False` AND
  `pytest.raises(MissingAssembly, match=<block message>)` on `grid.assemble()` (the block_array
  arm DEFERS — F3/R2). Tooth: monkeypatch `system_slices` to drop the System-B slice →
  round-trip reds. (This gate is NOT the leaf-type catcher — F2.)
- **G-c2.5 — forward block-`.H` reciprocity on the REAL grid (R1, the M4 row).** On the carrying
  sphere 2×2: `⟨grid.apply(ψ), x⟩_G = ⟨ψ, grid.H.apply(x)⟩_G` `< 1e-12` (grid.H is FREE via
  B.2a's `_AdjointOperator` + `CoupledSpace` member-wise metrics; measured 1.3e-16 on the toys —
  this is its FIRST real-curvilinear-member run). CONTROL holds + tooth **M-ADJ-metric** (skip
  the G⁺/G wrap in `_AdjointOperator.apply`, `operator.py:~1127`) → reds. A Euclidean block-`.H`
  reopens Mode-12/ERR-067. (The 1×1 non-carrying grid's `.H` is the existing FullFieldSpace
  composite adjoint — already Mode-12-closed; no new row.)
- **G-c2.6 — the dead-slot hazard-witness (R3, POSITIVE assert-defect).** Construct ψ_A with a
  LIVE ray = ψ_B; assert `grid.apply` DIFFERS from the fused by the seed-feed term (the
  double-count) — a POSITIVE assert-defect (like L16's out-of-scope-defect, NOT an xfail).
  Documents the transitional Pattern-4 hole; flips LOUD when B.2d's ray eviction makes a
  live-ray ψ_A unconstructable (the assert-defect must then be removed). Do NOT add a runtime
  guard forcing present-zero (calcifies a transient B.2d dissolves structurally).

---

## Mutation ledger (monkeypatch-only; NEVER `git checkout` an uncommitted file — L4)

| gate | tooth (in-process) | target `file:line` (⚠ best-estimate) | expected RED |
|---|---|---|---|
| G-c1.1 | drop `to_unified` bridge / swap interior↔boundary | NEW A_BB re-type `radial_characteristic.py:331/414` | container / `to_unified` value reds |
| G-c1.3 | emit bulk `AngularSourceSink` (drop FullField wrapper) | NEW A_AB re-type `radial_characteristic.py:734` | codomain != full_field_space / container reds |
| G-c2.1 (F2) | leave A_BB typed on the UNIFIED space in the builder | NEW `build_coupled_system` (B,B) wire | grid CONSTRUCTS (false green) but `grid.apply` CRASHES → runtime-apply gate reds |
| G-c2.2 | force System-B build on a slab (predicate always-true) | NEW `build_coupled_system` presence gate | block ctor raises `match=` its specific message |
| G-c2.3 (a) | swap A_AB↔A_BA placement | NEW builder block grid | O(1) centrepiece diff (asymmetric) |
| G-c2.3 (b) | drop `− B_b` from (B,B) | NEW builder (B,B) | ray rows red (reflective) / GREEN (vacuum — F1) |
| G-c2.4 | drop the System-B slice from `system_slices` | `numerics/coupled_system.py:483` (probe wrap) | round-trip reds |
| G-c2.5 | M-ADJ-metric (skip G⁺/G wrap) | `numerics/operator.py:~1127` (`_AdjointOperator.apply`) | reciprocity reds on sphere 2×2 |

---

## Re-point ledger (which existing tests MOVE, which STAY)

**MOVE (re-point to composite I/O — bit-id inheritance, `array_equal`):**
- `test_psi_half_coupling.py::TestA_BB_RadialBVP` (~9 foundation tests, ~15 ctor sites :799–:1200)
  — A_BB apply/solve/apply_transpose/solve_transpose/inverse probes → composite source/flux
  (build via `_ray_composite`; read `.to_unified()`/`.interior`). This IS G-c1.1/1.2 (F4).
- `test_ray_operator.py::TestA_BB_RadialBVP` (L1, 4 tests; the `_solve_…` helper :123) →
  composite source, read `.to_unified().cells`/`.interior`. Convergence orders UNCHANGED
  (bit-id engine).
- `test_psi_half_coupling.py::TestA_AB_SeedInjection` (~8 tests, :2299–:2600) — A_AB apply
  (seed→FullField, read `.interior`) / apply_transpose (FullField cotangent → composite, read
  `.to_unified()`). This IS G-c1.3 (F4; step-4 memo's "unaffected" was 4d.0-scoped).

**STAY (do NOT touch):**
- `tests/numerics/test_coupled_operator.py` (39, synthetic) — **0 ψ½/SN/build_coupled
  references** (grep-confirmed). The N-general machinery is CONSUMED by the builder, not
  changed. **NO changes.**
- `TestRegressionFloor` + `_dense`/`_blocks`/`_template`/`_loss` (FullField layout) — the fused
  `(L+C)` walk is UNTOUCHED at B.2c (step-4 refutation #7); NOT re-pointed until 4e.
- `TestA_BA_SchurFold`, `TestB_b_RayBoundary`, `TestCoupledLift` — already re-pointed at B.2b;
  B.2c does not re-touch them.
- `test_g_adjoint_reciprocity.py::_full_loss_case` — STAYS the reciprocity fixture (B_a-only is
  fine for reciprocity, NOTE R2). The centrepiece builds its OWN complete fused loss (F1).
- `test_inverse_adjoint_coherence.py` gates 1–3 — the swap-law `A.H.inverse()≡A.inverse().H` on
  `_loss(sn)` (the `(L+C)` op); the `coupled` fixture DEFERS to B.2d (R1).

---

## RULINGS on the open questions

- **R1 — A2a timing + E4 anchors: DEFER the inverse/swap-law + solve-anchors to B.2d; LAND the
  forward block-`.H` reciprocity at B.2c.** The A2a `coupled` fixture in
  `test_inverse_adjoint_coherence.py` (the swap law `A.H.inverse() ≡ A.inverse().H`) needs
  `grid.is_invertible`, which is `False` at B.2c (the block solve is step-5) — PREMATURE; land
  it at B.2d/step-5 with the block solve. But the FORWARD block-`.H` reciprocity (G-c2.5, NO
  inverse) IS cheap (grid.H free via B.2a) and load-bearing (the Mode-12/ERR-067 catcher's FIRST
  real-curvilinear-member run) — LAND it. E4 solve-anchors (φ=Q/Σ_t, k_inf) also defer (grid
  can't solve); the centrepiece's forward-equivalence to the already-anchored fused loss is the
  B.2c structural tie.
- **R2 — block_array/`assemble` arm DEFERS INDEFINITELY.** The ψ½ grid is a matvec/solve grid;
  its blocks don't emit sparse assembly (F3), so `grid.assemble()` raises `MissingAssembly`. The
  M2-on-real re-pin lands as `to_flat`↔`system_slices` coextensiveness + slice round-trip ONLY
  (G-c2.4), PLUS a gate asserting `is_assemblable is False` + `assemble()` raises (documents the
  unavailability so a future reader doesn't read it as a bug). The block_array arm re-enters only
  if a future SN walk-assembler wires `is_assemblable=True` on streaming — NOT in scope.
- **R3 — dead-slot: DOCUMENT the hazard + a POSITIVE assert-defect witness; do NOT guard.** A
  live-ray ψ_A silently double-counts the seed feed (LC(ψ_A) seeds the bulk on top of the
  explicit A_AB(ψ_B)). The centrepiece feeds ψ_A present-ZERO (well-posed) — Mode-12-BLIND to
  the double-count. Make it VISIBLE via G-c2.6 (a live-ray ψ_A → assert the double-count defect,
  positive, flips at B.2d). Document in the builder docstring. A runtime guard forcing
  present-zero would calcify a transient B.2d dissolves structurally (guard-at-source is B.2d's
  eviction, not B.2c's job).
- **R4 — A_BB solve-surface re-type scope = ALL FOUR action surfaces + involution-web re-point;
  `inverse()`/predicates untouched.** Re-type `apply`/`apply_transpose`/`solve`/`solve_transpose`
  (each bridges `to_unified`→engine→`from_unified`). `inverse()` returns `InverseOperator(self)`
  unchanged (it delegates to the re-typed solve/apply, so it speaks composite automatically).
  `is_invertible`/`is_adjointable` stay True. The involution-web gate (G-c1.2) re-points to
  composite I/O.

---

## Refutations (fire before the ink dries)

1. **F1: `_full_loss_case` omits B_b** — the centrepiece builds its OWN complete fused loss WITH
   B_b on a REFLECTIVE sphere (vacuum masks B_b). Do NOT reuse `_full_loss_case` as the
   centrepiece reference.
2. **F2: the block space-check is `==`-blind (unified ≡ composite by name+shape)** — "grid
   constructs" is necessary-NOT-sufficient (Mode-12); the sufficient catcher is the RUNTIME
   `grid.apply` (G-c2.1) + the centrepiece (G-c2.3). A_BB's re-type is a RUNTIME requirement.
3. **B.2c is NOT uniformly `array_equal`** — c1 re-types are (pure re-label, B.2b §0), but the
   c2 centrepiece bulk rows are principled-equiv `rtol=1e-11` (the intrinsic M-M block split,
   5.5e-16). A tight `array_equal` on `y_A.interior` falsely reds the honest split.
4. **P2 teeth are NET-NEW at the BUILDER** — the block ctors' seedless refusals exist, but the
   builder's presence-collapse (carrying→2×2 / non-carrying→1×1) has ZERO negative tests; write
   positive control + negative `match=`-specific.
5. **cyl/slab are the non-carrying CONTROL** — the builder degenerates to 1×1 `[[A_AA]]`; the
   offset-swap / centrepiece / reciprocity rows ride the carrying REFLECTIVE sphere-GL S4.
6. **the assemble arm is UNAVAILABLE (F3/R2)** — assert `is_assemblable False` + `assemble()`
   raises; the M2-on-real re-pin is `to_flat`↔`system_slices`, NOT `assemble≡probe`.
7. **the dead-slot is a documented hazard + positive witness (R3), NOT a guard** — flips at B.2d.

---

## Result contract

NEW file: `orpheus/sn/coupled_system.py` (`build_coupled_system`). NEW gate class:
`tests/sn/operators/test_psi_half_coupling.py::TestCoupledBuilder` (G-c2.1–c2.6). Extensions:
`TestA_BB_RadialBVP`/`TestA_AB_SeedInjection` in `test_psi_half_coupling.py` +
`TestA_BB_RadialBVP` in `test_ray_operator.py` (the c1 re-points = G-c1.1–c1.3). Every tooth
mutation-verified in-process under `-O`. End-to-end acceptance: `psi_half + ray_operator +
g_adjoint` green + full `tests/sn -m "not slow"` + `tests/numerics` (UNCHANGED, 39) + ratchet
`transport:1` + sphinx -W. **Load-bearing deliverables:** the **grid≡fused centrepiece
(G-c2.3, F1's complete fused loss, per-row tolerances)** + the **F2 runtime-apply proof
(G-c2.1)** + the **forward block-`.H` reciprocity (G-c2.5, the Mode-12 real-member catcher)**.
DEFER to B.2d: the A2a swap-law/inverse fixture, E4 solve-anchors, the dead-slot guard (R3
positive-witness now), the adapter retirement, the ray eviction.
