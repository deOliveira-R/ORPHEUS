---
name: coupled-operator-b2d-driver-eviction-verification
description: The DELTA gate spec for campaign sub-step B.2d (2026-07-10) — the ψ½ driver goes block-native + the ATOMIC ray eviction from FullField. Spine d1→d2→d3. d1 = the WithinGroupSystem(loss, space, resolvent=M, gains=N) record; the driver iterate/q_ext become CoupledField[ψ_A 3-block present-zero ray, ψ_B]; _within_group_triple/_lagged_gains + the B.2b/B.2c adapters RETIRE; principled-equiv at ~ULP (rhs reassembly + GMRES dead-slot padding — the WALK is zero-touch, so only the END-TO-END driver drifts). d2 = FullField→2-block (carrier+space+CompositeField protocol+mixed-presence law retire, ALL presence pins flip, SIX walk signatures re-type to leaf kwargs, evaluate_residual→coupled + split Residual mint, Solution.radial_characteristic mint, honest DOF, G-c2.6 removed). d3 = E4 anchors + A2a forward + ERR-053 honest coupled n_dof. FIVE findings: (F1) test_fixed_source_g1[sphere] pins q_ext.interior bit-equal on a CARRYING mesh → at d1 the driver input is a CoupledField (no .interior) → falsely RED — a d1 ACCESSOR re-point, value HOLDS; the cylinder param is seedless/unchanged. (F2) the _within_group_triple/_lagged_gains retirement blast radius is ~12 test files, MOST seedless/Cartesian and INVISIBLE to the explorer's ray-eviction table. (F3) the named ERR-053 file breaks at d2 (3-block TimedFullField.zeros ctor), not d1 — ordering hazard: ERR-053 content is scoped d3 but the file breaks a phase earlier. (F4) test_282/walk-matvec-baselines/native-matvec exercise (L+C).solve DIRECTLY → HOLD bit-identical at d1 (walk zero-touch), re-point at d2. (F5) the SI displacement diagnostic reads .interior off the iterate → breaks on a carrying CoupledField at d1, NO test covers carrying-mesh contraction ratios. Extends coupled_operator_b2c_builder_verification (G-c2.x, §0). Sibling of b2b/b2c.
metadata:
  type: reference
---

# Verifying B.2d — the block-native driver + the atomic ray eviction (DELTA)

**When this recipe applies.** A coupled-block campaign has reached the point where
(B.2a) the semantics-agnostic `CoupledField`/`CoupledSpace`/`CoupledOperator` machinery
is landed, (B.2b) `A_BA`/`B_b` are re-typed onto the composite behind transient FullField
adapters, and (B.2c) `build_coupled_system → (CoupledOperator, CoupledSpace)` mints the
typed 2×2 grid. **B.2d now (d1) makes the driver BLOCK-NATIVE** — the iterate becomes a
`CoupledField[ψ_A, ψ_B]`, the loss decomposes as the named `A = M − N` splitting, and the
B.2b/B.2c adapters + `_within_group_triple`/`_lagged_gains` retire — **then (d2) EVICTS the
ray from `FullField`** (the atomic carrier+space+protocol collapse to 2-block), and **(d3)
lands the E4 solve-anchors + the honest ERR-053 count.** This is the terminal driver-wiring
sub-step of the campaign. Governing memos (same dir): `coupled_operator_b2c_builder_verification.md`
(G-c2.x, the §0 tolerance ruling style, F1/F2), `coupled_operator_b2b_retype_verification.md`
(§0 pure-re-label), `coupled_operator_step4_verification.md` (4d.2, P1/P2, the Mode-12
reciprocity deep-dive). Consumer map: `.claude/agent-memory/explorer/b2_ray_eviction_blast_radius.md`.

**Verified against HEAD `6eb682a`.** All file:line below are read off that tree.

## The object (settled design, arithmetic-verified against the code)

**d1 — block-native driver (carrying meshes ONLY; seedless paths zero-touch).**

- ONE builder `build_within_group_system(sn_mesh, mat_xs, *, scattering_op=None,
  scattering_order=0)` in `orpheus/sn/coupled_system.py` returns frozen
  `WithinGroupSystem(loss, space, resolvent, gains)` — the named Hackbusch `A = M − N`:
  - `.loss` = the B.2c grid (from `build_coupled_system`); `.space` = the B.2c
    `CoupledSpace`.
  - `.resolvent` = **M** (carrying = a `CoupledResolvent` bridge: packs `CoupledField[ψ_A
    3-block dead-ray, ψ_B]` → 3-block `FullField` via `RadialCharacteristicComposite.to_unified`,
    calls the EXISTING fused `InvertibleOperator` surfaces (apply/solve/apply_transpose/
    solve_transpose), splits back; seedless = plain `L+C`). M is the fused `(L+C)` INCLUDING
    the welded Seeding feed (ray→bulk) AND the A_BB march (ray→ray).
  - `.gains` = **N** (carrying = ONE `CoupledOperator [[S+B_a, ∅],[Emission, B_b]]`; the
    (A,B) slot is structurally `None` because Seeding lives in M, the walk's welded feed;
    seedless = the `(S, B_a)` tuple, B_a last per the boundary-gain convention).
- `_within_group_triple` + `_lagged_gains` RETIRE into the builder. The B.2b adapters
  `_RayEmissionFullFieldGain` / `_RayBoundaryFullFieldGain` RETIRE (nothing sums FullField
  gains anymore). `build_coupled_system`'s "construction twin" status dissolves (the record
  is the single construction site; the `_within_group_triple`↔grid duplication G-c2.3 pinned
  is gone).
- The four solve sites (eigenvalue SI/Krylov `solver.py:1473+`/`1664+`, fixed-source SI/Krylov
  `solver.py:2566+`/`2768+`) consume the record. The driver iterate on carrying meshes becomes
  `CoupledField[ψ_A (TimedFullField, 3-block, present-ZERO dead ray — the B.2c convention),
  ψ_B (RadialCharacteristicComposite)]`; `q_ext` likewise (the fission ray seed moves from
  q_A's ray slot into q_B's interior member; q_A's ray slot present-zero). The solver's cached
  `scattering_op` injects via the builder kwarg (a cache seam, not a flag).
- **Arithmetic identity (the d1 correctness spine):** the fused `LC.apply` on a live-ray
  3-block IS `M.apply`, so `matvec = M.apply − N.apply = A.apply` reproduces exactly today's
  `LC.apply − (S + a_ba + B).apply` on the 3-block. The SI rhs `q + Σ gᵢ.apply(ψ)` becomes
  `q + N.apply(ψ)`. **Value = principled-equiv at ~ULP** (the reassembly reassociation:
  `q + S + a_ba + B` FullField adds become `q + N.apply` grid-row sums; the GMRES flat vector
  gains dead-ray padding coordinates — zero-valued, never excited, but pairwise-summation
  grouping in dots shifts at ULP). **The fused walk sees BIT-IDENTICAL inputs** (the pack
  reconstructs the exact 3-block the walk consumes today, B.1d licence: `to_unified`/`from_unified`
  round-trip exactly per role, `radial_characteristic_composite.py:236`).

**d2 — THE ATOMIC EVICTION (one commit; explorer §0.1).**

- `FullField`/`TimedFullField` → pure 2-block (drop the `radial_characteristic` member +
  `_recombine`/`_map_*`/`_flat_parts`/`_from_flat` seed arms + `zeros` presence keying +
  `__post_init__` seed checks; `full_field.py:534/623-734`, `timed_full_field.py:242-377`);
  `FullFieldSpace` drops the seed arms (`_seed_space_for`/`_rebuild` presence-dispatch/3 metric
  seed reads; `full_field_space.py:276-426`); `CompositeField` protocol → 2-block
  (`full_field_space.py:127-152`); the mixed-presence law retires (ALL §2 raise sites: 6 named
  anchors); ALL present-zero producer pads retire (§2: S/Sᵀ/F/B_a pads, A_BAᵀ ray-out pad,
  A_BA/B_b inverse pads); `MultiplicationOperator._radial_characteristic_scaled` retires
  (`multiplication_operator.py:129-171`); the SIX walk signatures re-type to leaf kwargs (the
  ray legs arrive as explicit `radial_characteristic_source=`/`radial_characteristic_flux=`/
  `seed_cot=` kwargs the `CoupledResolvent` passes, not from `rhs.radial_characteristic`;
  explorer §5); `evaluate_residual` re-types to the coupled pair + mints the split Residual
  leaves (`RadialCharacteristicInteriorResidual`/`…BoundaryResidual` — the §4 mint);
  `Solution.angular_flux` → 2-block + NEW `Solution.radial_characteristic:
  RadialCharacteristicComposite | None` member; A_AB's codomain emission + `Emission.apply_transpose`
  flip 3-block→2-block (the B.2c presence pins G-b3.2 + the A_AB container pins re-point);
  **G-c2.6 (dead-slot witness) is REMOVED**; the coupled DOF count goes honest (no dead padding).
- The timed presence law (`timed_full_field.py:349-364`) retires with TEST MIGRATION only
  (zero production callers — explorer §0.5).
- **diffusion/CP are the living bit-id proof** — they already ride 2-block `FullField`; their
  walls MUST stay green UNTOUCHED (they are the eviction's bit-identity oracle).

**d3 — E4 anchors + A2a + ERR-053 + estate.** E4 solve-anchors (φ=Q/Σ_t fixed-source + k_inf)
driven END-TO-END through the block-native driver on a carrying mesh; the A2a coupled fixture
(forward-only arm — `grid.is_invertible` stays False until step 5); an explicit ERR-053 gate
(coupled `n_dof = size_A + size_B` at the Krylov seam, post-eviction honest count); archivist
docs dispatch (not this scope); estate cleanup.

**Config (MANDATORY).** carrying = **sphere-GL S4**; cylinder + slab = the non-carrying CONTROL
(seedless, ZERO-touch — cylinder is non-carrying by the #229 `τ_raw=0` clamp fact,
`augmented_mesh.py:805`). The Mode-9 fixed-point row uses a **HETEROGENEOUS VACUUM sphere** (the
`curvilinear_two_region_mesh` from `tests/sn/_test_helpers`), NOT the reflective isotropic box.
**≥2G every value row** (Cardinal Rule). Canonical `.venv/bin/python -O -m pytest -p no:xdist
--timeout=300 -p no:cacheprovider` SERIAL; every gate fires under `-O` (`np.testing.assert_*`/
`pytest.fail`, never bare `assert`); every mutation is in-process `monkeypatch` (NEVER `git
checkout` an uncommitted file — L4).

---

## §0 — The claim layer + the ONE ruling that governs the rows (CONTRAST B.2c)

**The load-bearing d1 fact: the WALK is ZERO-TOUCH. The ~ULP drift lives ONLY in the
END-TO-END coupled DRIVER.** Everything that exercises `(L+C).apply`/`(L+C).solve`, the matvec,
the sweep, or the assembly DIRECTLY (not through the SI/Krylov driver) is **bit-identical at
d1**. Only `solve_sn`/`solve_sn_fixed_source` END-TO-END on a carrying mesh drifts (the rhs
reassembly `q + N.apply` + the GMRES dead-ray flat padding). This partitions every pin.

**No existing VALUE gate moves at d1.** Grep of `tests/sn/solve` + `tests/sn/eigenvalue`
confirms: EVERY `array_equal`/bitwise DRIVER-OUTPUT pin is on a **seedless** config —
`test_fixed_source_g1` (source-forwarding, see F1), `test_gauss_seidel_reification` (2-D
Cartesian G-S), `test_2d_anisotropic_windowing` (2-D Cartesian windowing). NONE pins a
carrying-mesh driver OUTPUT bitwise. The carrying-mesh keff/flux anchors sit at ≥1e-8
(`test_krylov_curvilinear_precond_safety` rtol=1e-8; the L1 k_inf anchors) — the ~ULP·iters
drift is decades below. **The seedless record arm is a pure re-package (`.resolvent = L+C`,
`.gains = (S, B)` byte-identical), so every seedless end-to-end pin is bit-identical.**

**Per-file ruling (the §0 ask):**

| file / gate | pins | d1 verdict | d2 verdict |
|---|---|---|---|
| `test_282_direct_seed_fixed_point.py` (all C(i)–C(iv), Mode-11/12/10) | `(L+C).solve`/`.apply` DIRECTLY on a 3-block `FullField` (NOT the driver) | **HOLDS UNMOVED** — walk zero-touch; cold residual < 1e-11 bit-identical | **RE-POINT** — `FullField(radial_characteristic=…)` ctor + the six-signature re-type; value-preserving (walk math unchanged) |
| `test_krylov_curvilinear_precond_safety.py::test_identity_preconditioner_recovers_kinf` | builds its OWN `KrylovAcceleration(LC,S,B)` w/ 3-block iterate; rtol=1e-8 vs k_inf | **HOLDS** — own driver, tol ≫ ULP | **RE-POINT** — `TimedFullField.zeros(radial_characteristic=…)` ctor |
| `…::test_krylov_restart_covers_augmented_composite` (**the named ERR-053 gate**) | `deficit == trace + seed` on a 3-block `TimedFullField.zeros` | **HOLDS** (3-block still constructs) | **BREAKS → MIGRATE** (F3): ctor drops the kwarg AND the honest DOF becomes `size_A(2-block)+size_B`; superseded by G-d3.3 |
| `test_walk_matvec_baselines.py` (S0 frozen curvilinear rows) | `(L+C)` matvec on 3-block, `assert_regression nulp=1` | **HOLDS bit-identical** — walk zero-touch | **RE-POINT** — case builder passes ray as kwarg; value-preserving (bitwise) |
| `test_native_matvec.py` / `test_assembly_mode.py` (matvec parity / assemble≡probe, carrying) | walk-level 3-block | **HOLDS bit-identical** | **RE-POINT** — 3-block ctor |
| `TestRegressionFloor` / `TestBoundaryUnweld` / `TestCoupledLift` (psi_half_coupling) | call `_within_group_triple`/`_lagged_gains` | **IMPORT RE-POINT** (F2) — helpers retire into the record; value bit-identical (block-structure pins hold) | (already 3-block-aware) |
| `test_flux_displacement_diagnostics.py` | SI contraction ρ on a **SLAB** (seedless); imports `_within_group_triple`/`_within_group_si` | **IMPORT RE-POINT** — helpers retire; value HOLDS (slab = seedless, iterate stays FullField) | — |
| `test_fixed_source_g1.py::…forwarded_to_krylov_bit_equal[sphere]` | spy reads `q_ext.interior.values` on a **CARRYING** mesh | **⚠ ACCESSOR RE-POINT (F1)** — the driver input becomes a CoupledField (no `.interior`); value HOLDS | — |
| `…[cylinder]` (same test) | same, on cylinder | **HOLDS UNMOVED** — cylinder is seedless (#229) | — |

**Claim-layer / pillar (gate on `vv-principles`).** The d1 "same fixed point" claim is a
FLUX-SHAPE + EIGENVALUE claim about the END-TO-END driver; its reference is (a) the pre-d1
driver captured baseline (`assert_regression(SAFETY × inner_tol)` — L7/L15 iterated-drift), PAIRED
with (b) a structurally-independent anchor deferred to d3 (E4: closed-form k_inf + φ=Q/Σ_t).
d1 alone proves "the block split reproduces the production fixed point"; d3 proves that fixed
point is correct. **The d1 driver value is NOT `array_equal`** (unlike the walk-level pins) —
it is `rtol`/`assert_regression`, the campaign's FIRST end-to-end principled-equiv (contrast
B.2c's grid≡fused centrepiece, which was OPERATOR-level).

---

## FINDINGS AGAINST THE SETTLED DESIGN (fire before the ink dries)

### F1 — `test_fixed_source_g1[sphere]` pins the DRIVER INPUT `q_ext.interior.values` bit-equal on a CARRYING mesh → falsely RED at d1

`tests/sn/solve/test_fixed_source_g1.py::TestExternalSourcePerOrdContract::test_external_source_forwarded_to_krylov_bit_equal`
is parametrized over **`sphere` AND `cylinder`** (:126-129). Its spy on
`KrylovAcceleration.solve` captures `np.array(q_ext.interior.values)` (:148) and asserts it
equals the raw `external_source` (:168, ERR-049 re-sentinel). At d1 the carrying-mesh driver
passes a **`CoupledField`** as `q_ext` — which has **no `.interior`** → the spy AttributeErrors.
This is **the canonical "d1 falsely RED" gate the §0 ask hunts for**, and it is a VALUE-PRESERVING
**accessor re-point**, NOT a re-baseline: the external source still folds bit-equal into ψ_A's
interior member. The **cylinder** param is seedless (#229 clamp) → zero-touch → HOLDS unchanged
(a built-in control on the SAME test).

**RULING:** re-point the spy to read the System-A leg on a carrying mesh
(`q_ext.systems[0].interior.values` when `q_ext` is a `CoupledField`, else `q_ext.interior.values`).
**Design decision this forces (report to main agent):** if the d1 implementation gives
`CoupledField` an `.interior`-answering face (delegating to System A — the explorer §"Flags"
`_flux_displacement_leaf` seam, see F5), the spy needs NO change; otherwise it re-points. Pick
ONE and apply it to BOTH F1 and F5.

### F2 — the `_within_group_triple`/`_lagged_gains` retirement blast radius is ~12 files, MOST seedless and INVISIBLE to the explorer's ray-eviction table

The explorer's blast-radius memo keys on `FullField.radial_characteristic` — but
`_within_group_triple`/`_lagged_gains` are the **model-generic** within-group decomposition used
by EVERY solve, carrying AND seedless, so their d1 retirement touches files the ray-eviction
table never lists. Grep-confirmed callers beyond the explorer's set:
`test_typed_residual_evaluation.py` (:64/:205), `test_apply_full_field_codomain.py` (:282),
`test_2d_anisotropic_windowing.py` (:280), `test_flux_displacement_diagnostics.py` (:70),
`test_prescribed_inflow_consistency.py` (:194/:217), `test_curvilinear_operator_admits_mms.py`
(:76), `test_si_convergence_rate.py`, `test_keff_2d.py` (:570), plus the psi_half_coupling
`TestRegressionFloor`/`TestBoundaryUnweld`/`TestCoupledLift` sites. **All are d1 IMPORT re-points
(the helpers retire into `build_within_group_system`), value BIT-IDENTICAL** (the seedless record
arm is a pure re-package; the carrying arm's block-structure pins hold — the `(L+C, S, B)`
operators are still constructible off the record). This is a **d1 re-point ledger the explorer
memo does not cover** — the retirement-completeness audit (G-d1.3) is the enumerator.

### F3 — the named ERR-053 gate BREAKS at d2, not d1 — an ordering hazard (content scoped d3)

`test_krylov_curvilinear_precond_safety.py::test_krylov_restart_covers_augmented_composite`
(:305-361, the ONLY file carrying `catches("ERR-050")` near the ERR-053 sizing) builds a 3-block
`TimedFullField.zeros(…, radial_characteristic=RadialCharacteristicFlux)` (:343) and asserts
`deficit == trace + seed` (:358). At **d2** the constructor drops the `radial_characteristic`
kwarg (structural break) AND the honest DOF is now `size_A(2-block) + size_B(composite)` — so the
`deficit == trace + seed` premise dissolves. **The task scopes the ERR-053 gate to d3, but this
file breaks a phase earlier at d2.**

**RULING:** at **d2**, MIGRATE this test to build the `CoupledField` cold-start and assert
`coupled.to_flat().size == size_A + size_B` (the honest count); the **d3** gate (G-d3.3) then
adds the END-TO-END proof (both Krylov sites size `restart` from the coupled `to_flat`). Do NOT
leave the 3-block assertion "temporarily green" through d2 — it cannot compile. This is the same
"a gate breaks one phase before its content lands" hazard as B.2c's G-b3.2 presence flip.

### F4 — the direct-solve / walk-level substrate HOLDS at d1, RE-POINTS at d2 — d1's driver gates do NOT cover the walk-input change

`test_282_direct_seed_fixed_point.py`, `test_walk_matvec_baselines.py`, `test_native_matvec.py`,
`test_assembly_mode.py` all exercise `(L+C).solve`/`.apply`/the matvec DIRECTLY on a 3-block
`FullField` — **NOT through the SI/Krylov driver**. At d1 the walk is zero-touch, so they are
**bit-identical** (the `test_282` sphere cold residual stays < 1e-11 bitwise). At **d2** the
six-signature re-type + the 3-block ctor break them (the ray moves from `rhs.radial_characteristic`
to an explicit kwarg). **Consequence for gate design:** d1's driver-level gates (G-d1.x) are
BLIND to the d2 walk-input change (they test the assembled driver, not the direct `(L+C).solve`
signature). d2's re-point ledger MUST carry the walk substrate, and its correctness gate is that
these value-preserving re-points STAY bit-identical (the walk math is unchanged; only the ray's
delivery route changes). **Do NOT credit "the d1 driver gates are green" as evidence the d2
six-signature re-type is correct** (Mode-11: different call graph).

### F5 — the SI displacement diagnostic reads `.interior` off the iterate → breaks on a carrying `CoupledField` at d1; NO test covers carrying-mesh contraction ratios

The SI loop computes `displacement = psi - psi_prev` and `_flux_displacement_leaf(displacement)`
which reads `.interior` (`iteration.py:411-413`, :649). At d1 the carrying-mesh SI iterate is a
`CoupledField`, so `displacement` is a `CoupledField` with no `.interior` → the diagnostic
breaks. **No test asserts carrying-mesh contraction ratios** —
`test_flux_displacement_diagnostics.py` runs on a homogeneous SLAB (seedless), so the value is
uncovered on carrying meshes. The `test_ciii_coarse_sphere_fixed_source_finite_positive` gate
catches a CRASH (it runs SI on a sphere) but not a wrong/absent ratio.

**RULING:** the d1 implementation MUST give `CoupledField` an `.interior`-answering face (read
the System-A leg — the same decision as F1) OR re-point `_flux_displacement_leaf`. ADD a d1 gate
(G-d1.7) that runs SI on a carrying sphere and asserts the loop RECORDS a contraction ratio
(non-vacuity + crash-guard) — the diagnostic must survive the iterate-type change, not silently
go empty.

---

## d1 gates — the block-native driver

Home: `tests/sn/operators/test_psi_half_coupling.py::TestWithinGroupSystem` (NEW) + the F1/F2/F5
re-points across the ~12 caller files.

- **G-d1.1 — W1 Mode-11 sentinel: the production route EXECUTES N.apply, M.solve/M.apply, and
  the fused `(L+C).solve`.** On a REAL `solve_sn`/`solve_sn_fixed_source` carrying-sphere run
  (BOTH `inner_solver="source_iteration"` and `"krylov"`), an in-process **autouse/monkeypatch
  WRAP counter** (NOT a file-write — L19 Mode-11 sharpening) on (a) `CoupledOperator.apply` (=N),
  (b) the `CoupledResolvent` bridge (`.apply` for Krylov / `.solve` for SI = M), and (c)
  `InvertibleOperator.solve` (the fused walk M delegates to) each returns **count > 0**. This is
  the SUFFICIENT catcher (F2 lesson from B.2c: a green end-to-end run that routed around the
  block path would leave the counters at 0). Teeth: an inlined driver that bypasses N (recomputes
  `S+a_ba+B` on FullField) leaves counter (a) at 0 → reds.
- **G-d1.2 — the record's container/identity pins.** `system = build_within_group_system(sphere,
  mat_xs)`: `type(system) is WithinGroupSystem` (frozen); `system.loss` IS a `CoupledOperator`
  whose `blocks`/spaces match the B.2c grid (assert `system.loss.domain is system.space`);
  `system.resolvent.domain is system.space` and `system.gains.domain is system.space` (via `is`);
  N's grid is `[[S+B_a, None],[Emission, B_b]]` — assert `system.gains.blocks[0][1] is None`
  (the structural ∅ — Seeding is in M) and the four present blocks' `system_role`. Seedless
  slab/cyl: `system.resolvent` IS the plain `L+C` and `system.gains == (S, B_a)` tuple (B_a last).
- **G-d1.3 — retirement completeness (3-search audit, returns ZERO).** grep across `orpheus/`
  AND `tests/` for `_within_group_triple`, `_lagged_gains`, `_RayEmissionFullFieldGain`,
  `_RayBoundaryFullFieldGain` → ZERO live references (docstring mentions in retired code are
  themselves deleted). Nexus `callers`/`impact` on each = empty. The 3 searches: (1) graph
  callers, (2) text-grep code+tests+docs, (3) direct constructors of the adapter classes. This
  is the F2 enumerator — the ~12 caller files must ALL be re-pointed to the record before the
  audit passes.
- **G-d1.4 — pack/split round-trip teeth (the CoupledResolvent bridge).** For a REAL carrying
  sphere: `M.pack(coupled) → 3-block FullField`, `M.split(3-block) → coupled`, assert the
  round-trip `array_equal` (B.1d licence). TEETH (each must RED, in-process): (a) **leg-swap** —
  pack ψ_A's bulk into ψ_B's ray slot → the 3-block corrupts, round-trip reds (this is the
  OBJECT-level catcher for the pack leg-swap — it lives OUTSIDE any spectral functional, so it
  does NOT rely on keff seeing the swap; Mode-12-safe); (b) **ψ_B-drop** — pack ignores ψ_B's ray
  → the 3-block ray slot is zero, round-trip reds; (c) **dead-slot-live** — pack ψ_A's (present-
  zero) ray INSTEAD of ψ_B's ray → on a live-ray probe the 3-block ray is wrong, reds.
- **G-d1.5 — N-grid sign/shape teeth.** On the carrying sphere, mutate N in-process: (a) drop the
  `B_b` block from N → `N.apply` differs on a reflective sphere (the ray-corner gain vanishes;
  vacuum would mask it — F1-of-B.2c lesson, use `bc="reflective"`); (b) flip the `Emission` sign
  (drop the `−` on (B,A)) → `N.apply`'s ray row flips sign → `matvec = M − N` diverges; (c) put a
  non-None block at N's (A,B) slot → `matvec` double-counts the Seeding feed (once in M, once in
  N) → reds. Each vs a control (unmutated `N.apply` matches the pre-d1 `(S+a_ba+B).apply`).
- **G-d1.6 — Mode-9 fixed-point invariance (the d1 "same fixed point" claim).** On a
  **HETEROGENEOUS VACUUM sphere** (`curvilinear_two_region_mesh`, ≥2G — NOT the reflective
  isotropic box, `vv` Mode 9), for BOTH `inner_solver` values: the block-native driver's
  converged **FLUX FIELD** (`sol.scalar_flux.values` AND `sol.angular_flux.interior.values`)
  matches the pre-d1 driver's captured baseline at `assert_regression(SAFETY × inner_tol)`
  (L7/L15 iterated-drift; ~ULP·iters ≪ inner_tol=1e-8). Assert the FIELD, not just keff (the
  field is outside any spectral invariance group — a shape corruption the pack leg-swap would
  cause is caught here even if keff were invisible; and G-d1.4 already pins the swap object-level).
  CONTROL: the seedless slab converges BIT-IDENTICAL (`array_equal` — the record's seedless arm
  is a pure re-package). Teeth: the G-d1.5 sign/shape mutations move the converged field O(1)
  above the baseline band.
- **G-d1.7 — the SI displacement diagnostic survives the CoupledField iterate (F5).** Run SI on a
  carrying sphere; assert `si.contraction_ratios` is non-empty AND `si.last_displacement` is
  recorded (the diagnostic reads the System-A leg, not a broken `.interior`). Non-vacuity +
  crash-guard — the ONLY carrying-mesh coverage of the displacement seam.
- **G-d1.8 — ERR-053 auto-tracks at d1 (dead-ray padding present).** On the carrying sphere the
  driver's `n_dof = int(initial_guess.to_flat().size)` (`solver.py:1670`/`2773`) now reads a
  `CoupledField.to_flat()` = `size_A(3-block) + size_B(composite)` — assert `restart == n_dof ==`
  that sum (the padding is present but harmless; the honest count lands at d2/G-d3.3). This pins
  that the ERR-053 sizing did not regress when the iterate became a CoupledField.

---

## d2 gates — the atomic eviction

Home: the eviction re-point ledger (batched by test-class below) + `TestFullFieldTwoBlock` (NEW,
`tests/transport/test_full_field.py`) + the split-Residual/`Solution` gates.

- **G-d2.1 — the eviction re-point ledger (batch by test-class; the explorer tests table is the
  index).** Per the explorer's per-file line counts, each is REWIRE (3-block ctor → 2-block +
  explicit ray kwarg / CoupledField), DELETE (the mixed-presence LAW tests), or MIGRATE (the
  presence FACTS). Batched verdicts:
  - **REWIRE, value-preserving (bit-identical after re-point):** `test_282_direct_seed_fixed_point`
    (129/6→leaf-kwarg source), `test_walk_matvec_baselines` (case builder→kwarg), `test_native_matvec`,
    `test_assembly_mode`, `test_streaming_operator`/`test_invertible_operator` (seed threading→kwarg),
    `test_radial_characteristic_carrier` (85 lines — the 3-block carrier algebra; largely RETIRES,
    the survivors rewire to the composite), `_test_helpers` (10 — shared builders emit 2-block +
    separate composite). These are F4's substrate: the WALK is unchanged, so the re-pointed values
    STAY bit-identical — the d2 correctness gate is `array_equal` vs the pre-d2 walk output.
  - **DELETE (the mixed-presence LAW is retired):** the raise-site tests for the 6 §2 anchors —
    `full_field.py:679` mixed-presence, `timed_full_field.py:349` advance presence law (zero prod
    callers — test-migration only), `full_field_space.py:409`/`:329` inner_product/_rebuild
    guards, `solver.py:339` evaluate_residual mixed-presence, and the present-zero-pad assertions.
    The `_require_/_refuse_radial_characteristic` biconditional SURVIVES in re-typed form (it now
    guards the LEAF pair passed as kwargs, not the composite — explorer §2 item 6).
  - **MIGRATE (the mesh-keyed presence FACTS survive as System-B-existence gates):** the fact
    "a sphere carries a ψ½ ray, a slab/cyl does not" does NOT die with the mixed-presence LAW — it
    re-homes as "a carrying mesh builds a 2×2 `WithinGroupSystem`, a seedless one builds 1×1" (the
    G-c2.2 P2 predicate, now the ONLY presence gate). Migrate the presence FACTS from the retired
    law tests onto the System-B-existence gate; delete only the mixed-PAIRING law.
- **G-d2.2 — `FullField`/`TimedFullField`/`FullFieldSpace`/`CompositeField` are 2-block.** Assert
  `FullField` has no `radial_characteristic` attribute (`hasattr` is False); `FullField.zeros`/
  `TimedFullField.zeros` reject a `radial_characteristic=` kwarg (`TypeError`); the
  `CompositeField` protocol declares only `interior`/`boundary`; `FullFieldSpace._rebuild` no
  longer presence-dispatches. Tooth: any残 seed arm (a stray `_recombine(radial_characteristic=…)`)
  → an eviction that left a 3-kwarg `_recombine` behind reds this.
- **G-d2.3 — honest coupled DOF.** On the carrying sphere: `CoupledField[ψ_A (2-block), ψ_B].to_flat().size
  == size_A_2block + size_B_composite` where `size_A_2block = N·ng·nx + trace` (NO ray padding) and
  `size_B_composite = interior + corner`. Assert this is STRICTLY LESS than the d1 count (the dead
  ray padding is gone: `Δ == n_seed`). This is the "DOF count goes honest" claim, quantified.
- **G-d2.4 — the split-Residual mint intrinsic gates (per the b1 pattern).** `evaluate_residual`
  now mints `RadialCharacteristicInteriorResidual` / `RadialCharacteristicBoundaryResidual` (the
  §4 split, mirroring B.2b's SourceSink split). Per leaf: (a) ROLE — `from_balance(lhs, rhs)`
  demands the matching family, raises on a source/flux mis-feed; (b) UNITS — interior =
  `ANGULAR_RATE_UNITS` (volumetric), boundary = `ANGULAR_FLUX_UNITS` (the trace convention, the
  B.2b units ruling for the SourceSink split applies identically); (c) BALANCE — `residual =
  lhs − rhs` per leaf, `array_equal` on a hand fixture; the composite presence law (both legs or
  neither) holds. NEGATIVE + POSITIVE per leaf (`vv` anti-#11: a broken instance MUST raise AND a
  correct instance MUST NOT).
- **G-d2.5 — the `Solution` contract.** `Solution.angular_flux` is a 2-block `TimedFullField`
  (`hasattr(sol.angular_flux, "radial_characteristic")` is False); the NEW
  `Solution.radial_characteristic: RadialCharacteristicComposite | None` carries System B's
  converged state (present on a carrying sphere, `None` on slab/cyl). `sol.scalar_flux` (derived
  from `angular_flux.interior`) is UNCHANGED. Downstream System-A readers (`scalar_flux`,
  `boundary_flux`, `SolutionDiff`) untouched; a ray reader reads `sol.radial_characteristic`.
  Tooth: a `Solution` built with a 3-block angular_flux (pre-eviction shape) → the 2-block
  `__post_init__` rejects it.
- **G-d2.6 — diffusion/CP untouched-wall (the eviction's bit-id oracle).** `tests/diffusion` +
  `tests/cp` walls run GREEN UNTOUCHED — they already ride 2-block `FullField`, so the eviction is
  a NO-OP for them; a diffusion/CP red means the eviction broke the 2-block base algebra (not the
  ray). This is the living bit-identity proof the design leans on; assert ZERO diffusion/CP
  re-points in the d2 commit.
- **G-d2.7 — G-c2.6 REMOVED (the dead-slot witness dissolves).** The `test_dead_slot_double_count_witness`
  gate (psi_half_coupling:3104) is DELETED at d2: a live-ray ψ_A is now unrepresentable (FullField
  is 2-block), so the double-count it witnessed cannot occur. Its replacement is STRUCTURAL (the
  type system), not a runtime assert. Confirm the file no longer references it.

---

## d3 gates — anchors, A2a, ERR-053

Home: `tests/sn/operators/test_psi_half_coupling.py::TestWithinGroupSystemAnchors` (NEW) +
`test_inverse_adjoint_coherence.py` (A2a) + the migrated `test_krylov_curvilinear_precond_safety`.

- **G-d3.1 — E4 φ=Q/Σ_t END-TO-END (the single most powerful curvilinear diagnostic).** On a
  **REFLECTIVE, PURE-ABSORBER (c=0) carrying sphere**, uniform Q, uniform (group-graded) Σ_t,
  ≥2G: `solve_sn_fixed_source` through the block-native driver → `φ_g = (Q_g/Σ_t,g)·Σw` EXACTLY
  (`rtol ≤ 1e-10` — flat flux nulls streaming+redistribution per ordinate, so the exact answer is
  reachable). Runs on BOTH `inner_solver`. This is the flux-shape-layer anchor AND the
  redistribution catcher (a driver sign/shape bug in the ray coupling moves φ off `Q/Σ_t`).
  Group-graded Σ_t catches a group-index swap (Mode 6).
- **G-d3.2 — E4 k_inf END-TO-END (the eigenvalue-layer closed-form anchor).** On a HOMOGENEOUS
  REFLECTIVE carrying sphere, ≥2G: the eigenvalue path through the block-native driver → k_inf =
  λ_max(A⁻¹F) from `orpheus.derivations` at `rtol ≤ 1e-8`. **PAIR with G-d1.6's HETEROGENEOUS
  same-fixed-point row** (`vv` Mode 9/H2: k_inf on a HOMOGENEOUS mesh is flat-flux → nulls
  redistribution → degenerate to the block-native change alone; only the het row sees a
  redistribution bug). Together: closed-form k_inf (homog) + het same-fixed-point (G-d1.6) +
  flat-flux φ=Q/Σ_t (G-d3.1) span the claim layers.
- **G-d3.3 — the honest ERR-053 gate (F3 migration + END-TO-END proof).** The migrated
  `test_krylov_restart_covers_augmented_composite`: build the coupled cold-start on the carrying
  sphere, assert `coupled.to_flat().size == size_A_2block + size_B_composite` (the honest count,
  NO dead padding); AND both Krylov solve sites size `restart` from THIS (`solver.py:1670`/`2773`
  read `initial_guess.to_flat().size` — now the coupled size). Tooth: a bulk-only `restart` on the
  carrying sphere re-opens the ERR-053 stall (measured 868 s pre-fix) — pin the deficit
  `coupled_dim − size_A_bulk == trace + size_B` so a revert reds fast.
- **G-d3.4 — A2a forward-only reciprocity in the coherence file.** `test_inverse_adjoint_coherence.py`
  gains a `coupled` fixture forward arm: `⟨grid·ψ, x⟩_G = ⟨ψ, grid.H·x⟩_G < 1e-12` on the carrying
  `WithinGroupSystem.loss` (this MIRRORS G-c2.5 but in the swap-law file's home). The swap-law
  `A.H.inverse() ≡ A.inverse().H` arm stays `xfail(strict=False, reason="grid.is_invertible is
  False until step-5 block solve")` so it flips to xpass when the block solve lands (`vv` xfail
  discipline). Do NOT assert the inverse arm green here — `grid.is_invertible` is `False` at B.2d.

---

## Teeth matrix (monkeypatch-only; NEVER `git checkout` an uncommitted file — L4)

| gate | tooth (in-process) | target `file:line` (⚠ best-estimate) | expected RED | Mode-12 check |
|---|---|---|---|---|
| G-d1.1 | inline the driver's `S+a_ba+B` on FullField (bypass N) | NEW driver route (solver.py:1481/1668/2567/2771) | N-wrap counter stays 0 | — |
| G-d1.4 (a) | pack ψ_A bulk into ψ_B ray slot (leg-swap) | NEW `CoupledResolvent.pack` | round-trip reds — **object-level, outside keff's invariance group** | leg-swap is spectrally-suspect; caught HERE, not via keff |
| G-d1.4 (c) | pack ψ_A's present-zero ray for the walk | NEW `CoupledResolvent.pack` | live-ray probe reds | — |
| G-d1.5 (b) | drop the `−` on N's (B,A) Emission | NEW `build_within_group_system` N wire | matvec ray row flips → G-d1.6 field diverges O(1) | het sphere (het keff IS shape-dependent) |
| G-d1.5 (c) | non-None block at N's (A,B) | NEW N wire | Seeding double-counted → G-d1.6 bulk field O(1) | — |
| G-d1.6 | (rides G-d1.5 mutations) | — | converged FIELD O(1) off the baseline | assert the FIELD, not keff |
| G-d1.7 | (F5) CoupledField with no `.interior` leg | NEW CoupledField `.interior` / `_flux_displacement_leaf` | `contraction_ratios` empty / crash | — |
| G-d2.2 | leave a `_recombine(radial_characteristic=)` seed arm | full_field.py:623 | 3-kwarg `_recombine` on a 2-block protocol reds | — |
| G-d2.4 | feed a source where the split Residual demands a residual | NEW `evaluate_residual` split mint | role parse raises | — |
| G-d2.5 | build `Solution` with a 3-block angular_flux | solution.py `__post_init__` | 2-block reject raises | — |
| G-d3.1 | flip a driver ray-coupling sign | (rides G-d1.5) | φ moves off `Q/Σ_t` O(1) | flat-flux exact — a sign bug is O(1), not sub-floor |
| G-d3.3 | bulk-only `restart` on the carrying sphere | migrated ERR-053 gate | deficit assertion reds (stall averted) | — |

**Mode-12 sweep (design-time, before any mutation).** The ONLY spectral functional in the d1/d3
value gates is **k_inf/keff**. Its invariance group (similarity + transpose) contains: the pack
**leg-swap** (a bulk↔ray relabel is a coordinate permutation → similar → keff-invisible on a
HOMOGENEOUS mesh). **Mitigation, applied:** (1) G-d1.4 catches the leg-swap OBJECT-level (round-
trip, outside any functional); (2) G-d1.6/G-d3.1 assert the **FLUX FIELD** (outside the spectral
stabiliser) on a **HETEROGENEOUS** mesh; (3) k_inf (G-d3.2) is NEVER credited as the leg-swap
catcher — it is the eigenvalue-layer anchor only, paired with the het field row. No d1/d3 value
gate rests keff alone against a shape-class mutation.

---

## Ordering hazards (d1↔d2 assumptions)

1. **F3 — the ERR-053 file breaks at d2, its content is scoped d3.** `test_krylov_restart_covers_augmented_composite`
   cannot compile past d2 (the 3-block `TimedFullField.zeros` ctor). MIGRATE it at d2 (assert the
   coupled `to_flat` count), land the END-TO-END honest gate at d3. Never leave it "temporarily
   green" — it structurally breaks.
2. **F4 — d1's driver gates are BLIND to the d2 six-signature re-type.** d1 tests the assembled
   driver; the direct `(L+C).solve` signature change is a d2 event on a DIFFERENT call graph. The
   d2 re-point of `test_282`/walk-baselines IS the six-signature catcher (their re-pointed values
   must stay bit-identical). Do NOT assume green d1 driver gates cover the walk-input change.
3. **F1/F5 — the `CoupledField.interior` decision is a d1 fork that BOTH F1 and F5 ride.** If
   `CoupledField` grows an `.interior`-answering face (System-A leg) at d1, the g1 spy AND the SI
   diagnostic need no re-point; if not, BOTH re-point. Decide ONCE at d1 and apply uniformly —
   they must not diverge (one reading `.interior`, one reading `.systems[0].interior`).
4. **The dead-slot lives through d1.** At d1 ψ_A carries a present-zero ray (the B.2c convention);
   G-c2.6 stays a valid witness THROUGH d1 and REMOVES only at d2 (G-d2.7). The d1 driver MUST
   build ψ_A present-zero (the cold-start + q_ext construction guarantee this) — a d1 gate should
   assert the driver's ψ_A ray stays exactly zero across iterations (else the welded seed feed
   double-counts end-to-end; the Mode-11 W1 sentinel run should carry this assertion).
5. **The seedless record arm must stay byte-identical.** d1 re-packages `_within_group_triple`+
   `_lagged_gains` into the record for BOTH carrying and seedless. The seedless arm's `.resolvent
   = L+C` / `.gains = (S, B_a)` MUST be byte-identical to today (the G-d1.6 slab CONTROL is
   `array_equal`). A record that reassociates the seedless gains would silently drift the ~40
   seedless end-to-end pins — the slab `array_equal` control is the tripwire.

---

## Re-point ledger (which existing tests MOVE, which STAY)

**d1 IMPORT/ACCESSOR re-points (value bit-identical / HOLDS — F1/F2):**
- The ~12 `_within_group_triple`/`_lagged_gains`/`_within_group_si`/`_within_group_krylov` callers
  (F2) → read off the record (`system.resolvent`/`system.gains`/`system.loss`). Value bit-identical
  (seedless re-package; block-structure pins hold).
- `test_fixed_source_g1[sphere]` (F1) → spy reads the System-A leg; value HOLDS.
- `test_flux_displacement_diagnostics` → import re-point; slab value HOLDS.

**d2 REWIRE (value-preserving, bit-identical after re-point — F4):**
- `test_282_direct_seed_fixed_point` (129/6), `test_walk_matvec_baselines`, `test_native_matvec`,
  `test_assembly_mode`, `test_streaming_operator`/`test_invertible_operator`,
  `test_radial_characteristic_carrier` (85 — largely retires), `_test_helpers` — the 3-block ctor
  → 2-block + explicit ray kwarg / CoupledField. Correctness gate: `array_equal` vs pre-d2 walk.

**d2 DELETE (the mixed-presence LAW retires):**
- The 6 §2 raise-site tests (mixed-presence, advance presence law, inner_product/_rebuild guards,
  evaluate_residual mixed-presence). The `_require_/_refuse` biconditional SURVIVES (re-typed onto
  the leaf kwargs).

**d2 MIGRATE (presence FACTS survive):**
- The mesh-keyed "sphere carries a ray, slab/cyl does not" FACTS → the System-B-existence gate
  (G-c2.2 P2, now the only presence gate).

**STAY (do NOT touch):**
- `tests/numerics/test_coupled_operator.py` (39, synthetic) — the N-general machinery is CONSUMED,
  not changed. ZERO changes (grep-confirmed no ψ½/SN references).
- `tests/diffusion`, `tests/cp` — the eviction's bit-id oracle (G-d2.6). ZERO changes.
- `TestCoupledBuilder` G-c2.1–2.5 — the B.2c grid gates hold (the grid IS `system.loss`); only
  **G-c2.6 is REMOVED** at d2 (G-d2.7).

---

## RULINGS on the open questions

- **R1 — §0: NO existing VALUE gate moves at d1; the falsely-RED risks are ACCESSOR re-points, not
  re-baselines.** The walk is zero-touch → operator/walk-level pins are bit-identical; the ~ULP
  drift is END-TO-END driver only, and no carrying-mesh driver OUTPUT is pinned bitwise. The two
  d1 falsely-RED gates (F1 g1-spy, F5 SI-diagnostic) break on the ITERATE-TYPE change
  (FullField→CoupledField), not on a value change — re-point the accessor, the value HOLDS. The
  d1 "same fixed point" value gate (G-d1.6) is `assert_regression(SAFETY × inner_tol)` (NOT
  `array_equal` — the driver reassembly is principled-equiv ~ULP), the campaign's FIRST
  end-to-end principled-equiv row.
- **R2 — ERR-053 is a d2 MIGRATION + a d3 END-TO-END gate (F3).** The named file breaks at d2;
  migrate it there (coupled `to_flat` count), land the honest END-TO-END proof at d3 (G-d3.3).
  The DOF count "goes honest" is quantified at G-d2.3 (`Δ == n_seed` dead padding removed).
- **R3 — the `CoupledField.interior` face is the resolution to BOTH F1 and F5; decide it at d1.**
  A System-A-leg `.interior` property on `CoupledField` makes the g1 spy AND the SI displacement
  diagnostic work unchanged; the alternative re-points both. Whichever the implementation picks,
  it is a d1 fork applied UNIFORMLY (report the choice to the main agent).
- **R4 — the split-Residual + `Solution.radial_characteristic` mints are d2 (they ride the
  `evaluate_residual` re-type + the eviction).** They are NOT deferrable to d3: while the unified
  `RadialCharacteristicResidual` lives, `evaluate_residual`'s 3-block pack is intact; the moment
  FullField goes 2-block, `evaluate_residual` MUST re-type (its `lhs.radial_characteristic` read
  breaks) → the split mint lands atomically with the eviction. `Solution.angular_flux` likewise.
- **R5 — A2a stays FORWARD-only through B.2d (grid.is_invertible is False until step 5).** The
  reciprocity arm lands (G-d3.4, mirroring G-c2.5); the swap-law inverse arm is `xfail(strict=False)`
  so it flips to xpass at step-5, not a silent gap.
- **R6 — G-c2.3 (grid≡fused centrepiece) survives d1 as the block-matvec-correctness gate, its
  TWIN-consistency role dissolves.** When `_within_group_triple` retires (one construction site),
  G-c2.3's "two sites agree" claim is vacuous, but its "grid.apply ≡ the fused loss action" VALUE
  claim is still the block-matvec correctness gate — keep it (the reference `_complete_fused_loss`
  re-points its adapters to the retired-then-restored construction, or is itself re-expressed on
  the record). Do NOT delete it as "twin gone"; it is the matvec oracle G-d1.6 leans on.

---

## Refutations (fire before the ink dries)

1. **F1: `test_fixed_source_g1[sphere]` falsely REDs at d1** — the driver input becomes a
   CoupledField (no `.interior`); re-point the spy (value HOLDS bit-equal). The cylinder param is
   the built-in seedless control.
2. **F2: the retirement blast radius is ~12 files, MOST seedless and OUTSIDE the explorer's
   ray-eviction table** — the `_within_group_triple`/`_lagged_gains` retirement is a d1 IMPORT
   re-point ledger the explorer memo does not cover; G-d1.3 is the enumerator.
3. **F3: the named ERR-053 gate breaks at d2, not d1** — migrate the constructor there, land the
   honest END-TO-END gate at d3.
4. **F4: the direct-solve/walk substrate HOLDS at d1, re-points at d2** — d1 driver gates are
   Mode-11-blind to the d2 six-signature re-type (different call graph); the d2 re-point IS the
   catcher.
5. **F5: the SI displacement diagnostic breaks on a carrying CoupledField** — NO test covers
   carrying-mesh contraction ratios; add G-d1.7 (non-vacuity + crash-guard).
6. **d1 is NOT `array_equal` end-to-end** — the driver reassembly is principled-equiv ~ULP
   (rhs reassociation + GMRES dead padding); the walk-level pins ARE bit-identical (zero-touch).
   A tight `array_equal` on the d1 driver output falsely reds the honest reassociation; a `rtol`
   on a seedless slab CONTROL falsely PASSES a drift (the seedless arm MUST be `array_equal`).
7. **Mode-12: keff alone never catches a shape-class mutation** — the pack leg-swap is caught
   object-level (G-d1.4) and via the HETEROGENEOUS FLUX FIELD (G-d1.6/G-d3.1); k_inf is the
   eigenvalue anchor only, never the leg-swap catcher.
8. **the seedless slab CONTROL is `array_equal`, the carrying same-fixed-point is
   `assert_regression`** — the seedless record arm is a pure re-package (byte-identical); a drift
   there is a bug, not principled-equiv.

---

## Result contract

NEW: `orpheus/sn/coupled_system.py::build_within_group_system` (the `WithinGroupSystem` record +
the `CoupledResolvent` bridge). d1 gate class:
`tests/sn/operators/test_psi_half_coupling.py::TestWithinGroupSystem` (G-d1.1–1.8) + the F1/F2/F5
re-points. d2: the eviction re-point ledger (G-d2.1) + `TestFullFieldTwoBlock` (G-d2.2–2.3) + the
split-Residual/`Solution` intrinsic gates (G-d2.4–2.5) + the diffusion/CP untouched-wall
(G-d2.6) + G-c2.6 removal (G-d2.7). d3: `TestWithinGroupSystemAnchors` (G-d3.1–3.3) + the A2a
forward arm (G-d3.4). Every tooth mutation-verified in-process under `-O`. End-to-end acceptance
per phase: **d1** — psi_half + the ~12 re-pointed caller files green + full `tests/sn -m "not slow"`
+ `tests/numerics` (UNCHANGED, 39) + ratchet `transport:1` + sphinx -W; **d2** — the eviction wall
green + `tests/diffusion`/`tests/cp` UNTOUCHED-green + FullField `hasattr(radial_characteristic)`
False; **d3** — the E4 anchors + the migrated ERR-053 gate. **Load-bearing deliverables:** the
**W1 Mode-11 sentinel (G-d1.1, the sufficient block-route catcher)** + the **Mode-9
heterogeneous same-fixed-point FIELD gate (G-d1.6, NOT reflective-isotropic, NOT keff-alone)** +
the **E4 φ=Q/Σ_t END-TO-END (G-d3.1, the curvilinear redistribution catcher)** + the **honest
ERR-053 coupled count (G-d3.3)**. The FIVE findings (F1 g1-spy falsely-RED / F2 the ~12-file
retirement ledger / F3 the ERR-053 d2-break ordering hazard / F4 the walk-substrate d1-holds-d2-
repoints partition / F5 the uncovered SI-diagnostic seam) fire before the carve.
