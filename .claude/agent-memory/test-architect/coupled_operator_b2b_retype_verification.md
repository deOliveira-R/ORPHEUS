---
name: coupled-operator-b2b-retype-verification
description: The DELTA gate spec for campaign sub-step B.2b (2026-07-10) — re-type A_BA (codomain) + B_b (domain/codomain) onto RadialCharacteristicComposite behind thin transient FullField-gain adapters (production BIT-IDENTICAL), minting the split-locus SourceSink pair the block codomains need. Governs three sub-commits (b1 leaves+bridge / b2 member space / b3 operator re-type). THE ruling: B.2b is a PURE RE-LABELING — every value row is array_equal (bit-id inheritance), NO principled-equiv row (rtol/nulp in a B.2b gate is a RED FLAG). Extends the step-4 memo's P3 (value-unmoved) + re-point ledger. Sibling of coupled_operator_step4_verification (governing) + coupled_operator_b1_split_verification (B.1d leaves/composite).
metadata:
  type: reference
---

⛔ **Displacement/FluxRole rows OVERTURNED (campaign 1 CS3, 2026-08-19).**
The role set below is now {Flux, SourceSink} — the Displacement role and the
`FluxRole` mixin retired with the cone carve (flux lives in V), so the
Displacement rows of the bridge/role tables and the "inherit `FluxRole`"
mutation are historical. The SourceSink minting, adapters, and bit-identity
rulings stand.

# Verifying B.2b — the A_BA / B_b re-type onto System B's composite (DELTA)

**When this recipe applies.** A coupled-block campaign has (B.1) posed System B as
an independent `RadialCharacteristicComposite = Composite[interior ⊕ boundary]` and
(B.2a) landed the semantics-agnostic `CoupledField`/`CoupledSpace`/`CoupledOperator`
machinery. B.2b now **re-types the two System-B-touching SN operators onto that
composite** while keeping production BIT-IDENTICAL through thin transient
FullField→FullField gain adapters. It is the RE-POINT half of the step-4 memo's P3;
this delta memo is the gate spec for its three sub-commits. Governing memo:
`coupled_operator_step4_verification.md` (P3 + carry-forward/re-point ledger).
Prior leaf/composite gates: `coupled_operator_b1_split_verification.md`.

**The object.** Three sub-commits, additive:
- **b1** — mint the split-locus SOURCE-SINK role pair
  `RadialCharacteristicInteriorSourceSink` + `RadialCharacteristicBoundarySourceSink`
  (over the `_bases.py:1242/1304` field bases, `ANGULAR_RATE_UNITS`, mirroring the
  unified `RadialCharacteristicSourceSink`); re-bind `RadialCharacteristicComposite`'s
  **static** type params from the flux leaves to the FIELD BASES (role-erased — the
  FullField `Composite[BulkField, BoundaryField]` precedent; the runtime `__post_init__`
  guard ALREADY admits bases); make the bridge `from_unified`/`to_unified` ROLE-PRESERVING
  via an exact-class leaf table (Flux↔(InteriorFlux,BoundaryFlux); SourceSink↔the new
  source pair; Displacement↔the split displacement pair).
- **b2** — System B's member SPACE: `FullFieldSpace.from_blocks` gains a `name` param
  (System B named `"radial_characteristic"`, built from `interior=ray_interior_space,
  trace=ray_boundary_space`); `FullFieldSpace._rebuild` (full_field_space.py:300-344)
  gains PRESENCE-DISPATCH (pass the `radial_characteristic` kwarg to `_recombine` ONLY
  when the carrier exposes the slot — today it unconditionally reads `x.radial_characteristic`
  AND passes the kwarg, which `AttributeError`s / `TypeError`s on a pure 2-block
  `RadialCharacteristicComposite`); new cached `SNMesh.radial_characteristic_composite_space`.
- **b3** — re-type `A_BA` = `RadialCharacteristicEmission` (codomain → the ray composite;
  `apply` returns a `RadialCharacteristicComposite` with SOURCE members, no FullField
  present-zero padding; `apply_transpose` takes the composite cotangent, returns FullField
  with `radial_characteristic=None`, the forward-looking 2-block System-A shape); re-type
  `B_b` = `RadialCharacteristicBoundaryOperator` (domain=codomain → the ray composite;
  reads the boundary member's flux corners, emits boundary-member source corners). The
  fold/reflect engines stay UNIFIED internally behind the bridge (the demote ruling —
  full unified-leaf retirement is Phase C/4e). Production stays BIT-IDENTICAL via
  solver-private FullField→FullField adapters at the two gain seams (`_lagged_gains`
  wraps A_BA; the `B_a + B_b` sum in `_within_group_triple` wraps B_b).

**Carrying member = sphere-GL S4 ONLY (R12a).** cylinder-LS + slab = non-carrying CONTROL
— System B's leaves/composite/member-space are UNCONSTRUCTABLE there (the presence-raise IS
the control gate, not a degenerate value). **≥2G every value row** (Cardinal Rule; 1G k is
flux-shape-independent). Canonical invocation
`.venv/bin/python -O -m pytest -p no:xdist --timeout=300 -p no:cacheprovider` SERIAL; every
gate fires under `-O` (`np.testing.assert_*` / `pytest.fail`, never a bare `assert`); every
mutation is in-process monkeypatch (NEVER `git checkout` an uncommitted file — L4).

## §0 — The claim layer, the pillar, and the ONE ruling that governs all rows

**Claim layer (gate on `vv-principles`): B.2b is a PURE RE-LABELING.** The re-type moves
NO arithmetic — A_BA computes the identical `Fold ∘ K_iso ∘ integrate`, then bridges the
unified ray leaf into the composite (a pure gather/scatter, `from_unified`); B_b applies the
identical specular corner swap on the boundary member. So every value row is a **bit-identity
claim** (equivalence bar = `array_equal`, verified by INHERITANCE from the untouched kernel).
The eigenvalue/flux-shape claims (E4 `φ=Q/Σ_t`, `k_inf`) DO NOT MOVE and are NOT re-proven
per gate — they are the end-to-end wall the existing `tests/sn` suite owns.

**THE RULING — `array_equal` everywhere; a `rtol`/`nulp` in a B.2b gate is a RED FLAG.**
Because the re-type touches no reduction tree, there is NO principled-equivalence row in
B.2b (contrast the un-weld/extract steps where COO→CSR / iterated-drift earns `rtol=1e-11`).
The FIRST principled-equiv drift enters at B.2d/4e. If any B.2b production-value gate needs a
tolerance to pass, the re-type moved arithmetic it should not have — investigate, do NOT
relax (coding-elegance anti-#17 / vv §Bit-identity).

**Pillar / structural independence.** The bit-id inheritance reference for the A_BA value is
`_ba_oldloop_reference` (the documented fold loop — the retained fuller-view oracle sharing
the fold kernel; representation differs: unified leaf vs composite), and for B_b the
hand-spelled specular swap (`corner(level,−1) = seed.corner(level,+1)`). Both survive the
carve (they are NOT on the deleted old-operator path), so they ARE the frozen pre-carve
reference in spirit (L19). Belt-and-suspenders option: a root-conftest `--capture-baseline`
snapshot captured BEFORE b3. The lower-layer structural anchor (φ=Q/Σ_t, k_inf) is already
discharged pre-B.2b; B.2b sits ABOVE it as re-labeling — do NOT re-anchor it here.

---

## b1 — the SourceSink pair + role-preserving bridge + re-bound composite

### Net-new gates (extend `tests/sn/mesh/test_radial_characteristic_split_leaves.py` [leaves] + `tests/transport/test_radial_characteristic_composite.py` [bridge/composite])

- **G-b1.1 — source-leaf intrinsic laws (cross-role sum rejection + closed source algebra).**
  For BOTH `RadialCharacteristicInteriorSourceSink` and `…BoundarySourceSink`:
  (i) POSITIVE — `source + source` → a source of the SAME class, `array_equal(values, a+b)`
  (SourceSinks are plain vector-space, no role mixin — source sums closed);
  (ii) NEGATIVE — `source + <matching-locus Flux>` → `TypeError` (`match=` the affine/role
  message) — a cross-role mix is unrepresentable. `pytest.raises` with `match=` (a bare
  `raises` false-greens on a downstream shape crash — anti-#11).
  - **bit-id ruling:** the `+` VALUE is `array_equal` (exact).
  - **Mode-12 lens:** the `type(a+b) is …SourceSink` assertion is blind to the buffer VALUE
    (invariant under any permutation of `.values`); PAIR it with `array_equal(values)`. Each
    is blind to the other's error class.
  - **Teeth (in-process):** monkeypatch the SourceSink to inherit `FluxRole` (so `source+source`
    would mint a displacement / raise) → the POSITIVE control reds; the cross-role negative is
    self-toothed by its `match=`.

- **G-b1.2 — source-leaf presence-raise on a non-carrying mesh.** `…InteriorSourceSink.zeros_on(slab)`
  → `ValueError(match="carries no radial_characteristic_interior_space")`; same for the boundary
  sibling on `…radial_characteristic_boundary_space`. POSITIVE control: `.zeros_on(sphere)` succeeds
  (nonzero-shaped). Mirrors the existing flux-leaf presence gate (`_face_space_of` raises).
  - **bit-id ruling:** structural (constructability), not a value.
  - **Mode-12 lens:** presence is a structural predicate — no functional blindness; the `match=`
    pins the SPECIFIC absent-space message (a bare `raises` false-greens on any ValueError).
  - **Teeth:** the negative IS the tooth (paired with its positive control — anti-#11).

- **G-b1.3 — THE b1 HIGHEST-VALUE ROW: bridge role-preservation, per role, bitwise.** For each
  role `u ∈ {Flux, SourceSink, Displacement}`, build a unified `u` (random, nonzero), then
  `c = RadialCharacteristicComposite.from_unified(u)`:
  (i) ROLE — `type(c.interior) is <matching InteriorRole>` AND `type(c.boundary) is <matching
  BoundaryRole>` (Flux→InteriorFlux/BoundaryFlux; SourceSink→the new pair; Displacement→the
  split displacement pair);
  (ii) VALUES — `array_equal(c.to_unified().values, u.values)` (round-trip bitwise, cells→interior
  + corner→boundary).
  - **bit-id ruling:** `array_equal` (the bridge is a pure gather/scatter).
  - **Mode-12 lens (THE split):** the value round-trip is blind to the leaf CLASS — the CURRENT
    hardcoded bridge (always mints Flux leaves) would PASS `array_equal` for ALL three roles yet
    be WRONG for SourceSink/Displacement. And a value-scrambling bridge would pass the class check
    yet fail `array_equal`. BOTH assertions mandatory per role.
  - **Teeth + necessity (§0.6 asymmetry):** monkeypatch the bridge's role-table so SourceSink→Flux
    (the pre-b1 hardcode) → the ROLE assertion reds for the SourceSink AND Displacement rows while
    STAYING GREEN for the Flux row — the asymmetry IS the config-blindness evidence (a Flux-only
    bridge test is structurally blind to the role-preservation bug). Second tooth: a bridge that
    threads cells but drops the corner → the `to_unified` boundary round-trip reds.

- **G-b1.4 — the re-bound composite ADMITS source members.** `RadialCharacteristicComposite(
  interior=…InteriorSourceSink.zeros_on(sn), boundary=…BoundarySourceSink.zeros_on(sn))`
  constructs WITHOUT raising (the runtime `__post_init__` guard on the field BASE admits the
  source role). POSITIVE control: a flux-member composite still constructs. NEGATIVE: a foreign
  interior (`AngularFlux`) still raises `TypeError(match="RadialCharacteristicInteriorField")`.
  - **Static half (pyright, NET-NEW):** `assert_type(composite.interior, RadialCharacteristicInteriorField)`
    and a `# type: ignore`-free source-member construction — the static-param re-bind from the
    flux leaf to the field base is what makes it type-check. This is a PYRIGHT gate: reverting the
    type param to `RadialCharacteristicInteriorFlux` reveals `reportArgumentType` on the
    source-member construction (the RUNTIME guard already admits it — the STATIC gate is the catcher).
  - **bit-id ruling:** structural.
  - **Mode-12 lens:** construction-success is structural; the pyright `assert_type` is the static
    contract. Runtime + static are complementary — a runtime-only gate is blind to the static-param
    regression.

### Re-points (b1): NONE — b1 is purely additive (new leaves, a wider composite type param, a role-preserving bridge). The existing B.1d bridge round-trips (`TestSplitFidelityBridge`) stay GREEN and byte-identical for the Flux role (the additive-migration contract); G-b1.3 EXTENDS them to the two new roles.

---

## b2 — System B's member space (FullFieldSpace reuse + presence-dispatch)

### Net-new gates (new file `tests/sn/mesh/test_radial_characteristic_composite_space.py`, mirroring the split-space tests)

- **G-b2.1 — member-space identity + metric-trio ≡ direct split-space application on the REAL
  composite.** `space = sn.radial_characteristic_composite_space`:
  (i) IDENTITY — `isinstance(space, FullFieldSpace)`, `space.name == "radial_characteristic"`,
  `space.interior_space is sn.radial_characteristic_interior_space`, `space.trace_space is
  sn.radial_characteristic_boundary_space`, `space.radial_characteristic_space is None` (pure
  2-block);
  (ii) METRIC-TRIO on a REAL `RadialCharacteristicComposite` `c` (random flux members):
  `space.apply_metric(c)` ≡ the direct split-space application (interior member scaled by
  `interior_space.apply_metric`, boundary by `trace_space.apply_metric`) — `array_equal` PER BLOCK
  (object-level); `space.apply_inverse_metric(c)` ≡ direct — `array_equal`; `space.inner_product(c, c')`
  == `interior_space.inner_product(...) + trace_space.inner_product(...)` (scalar).
  - **bit-id ruling:** `array_equal` (the metric arithmetic is untouched — only presence-dispatched).
  - **Mode-12 lens:** `space.__eq__` is `(name, shape)` only (`interior_space`/`trace_space` are
    `compare=False`) → `==` is BLIND to the member spaces; assert `interior_space is` / `trace_space is`
    the split spaces DIRECTLY. And `inner_product` is a SCALAR functional (blind to any permutation
    inside the metric's stabiliser) — PAIR it with `apply_metric`/`apply_inverse_metric` `array_equal`
    (which pin the per-block placement, the OBJECT).
  - **Teeth:** (a) PRESENCE-DISPATCH INVERSION — monkeypatch `_rebuild` to unconditionally read
    `x.radial_characteristic` + pass the kwarg (the pre-b2 behavior) → `apply_metric(real_composite)`
    raises (`AttributeError`/`TypeError`: the 2-block `_recombine` rejects the kwarg) — the metric-trio
    reds. (b) FullField BYTE-IDENTITY — `apply_metric(fullfield_with_seed)` UNCHANGED (`array_equal` vs
    the pre-b2 value on a seed-carrying FullField); tooth: a presence-dispatch that dropped the seed
    kwarg for FullField too → the FullField ψ½-block metric vanishes → reds. **BOTH the FullField-seed
    path AND the 2-block path must be exercised** — the presence-dispatch is exactly the branch that
    discriminates them, so a suite testing only one is blind to the wrong branch of the dispatch (L11).

- **G-b2.2 — the coextensiveness pin `to_flat().size == prod(space.shape)` on the REAL ray composite.**
  `RadialCharacteristicComposite.from_mesh(sphere).to_flat().size == int(np.prod(
  sn.radial_characteristic_composite_space.shape))`. (The B.2a follow-up reminder made concrete:
  the toys had this by contract; on the REAL multi-axis member the three offset spellings —
  `to_flat`, `from_blocks` shape sum, the member metric layouts — must AGREE.)
  - **bit-id ruling:** structural.
  - **Mode-12 lens:** a size EQUALITY is blind to layout ORDER (two layouts of equal total pass);
    companions that pin order already exist — the B.1d `from_flat`/`to_flat` round-trip + G-b2.1's
    per-block `apply_metric` `array_equal`.
  - **Teeth:** monkeypatch `from_blocks` to drop the boundary block from the shape sum → the
    coextensiveness reds.

### Re-points (b2): NONE — b2 adds a `name` param (default preserves `"full_field"` byte-identity for System A) + a new mesh prop + the `_rebuild` presence-dispatch. The `_rebuild` change is the ONLY behavior touch, and its FullField-neutrality is G-b2.1 tooth (b). The existing FullFieldSpace suite (`test_full_field_space.py`, `test_g_adjoint_reciprocity.py`, `test_starting_direction_metric.py::test_shipped_metric_block_values` atol=0.0) MUST stay GREEN + UNMODIFIED (the FullField metric path is byte-identical — the presence-dispatch only ADDS the 2-block branch).

---

## b3 — the operator re-type (A_BA codomain, B_b domain=codomain) + adapters

### Net-new gates (extend `tests/sn/operators/test_psi_half_coupling.py`)

- **G-b3.1 — A_BA block apply ≡ old embedded value (P3, bit-id INHERITANCE).** `out = A_BA.apply(psi)`
  (sphere, ≥2G, nonzero `psi`): (i) VALUE — `array_equal(out.to_unified().values,
  _ba_oldloop_reference(emission, sn))`; (ii) CONTAINER — `type(out) is RadialCharacteristicComposite`,
  `type(out.interior) is RadialCharacteristicInteriorSourceSink`, `type(out.boundary) is
  …BoundarySourceSink`; (iii) the OLD "bulk/boundary present-zero" is now STRUCTURAL — A_BA's codomain
  has NO bulk slot, so "A_BA writes into the bulk (double-count with S)" is UNSPELLABLE (Pattern-4).
  - **bit-id ruling:** BIT-IDENTICAL — the fold arithmetic is textually unchanged; A_BA computes the
    unified ray then `from_unified` (pure gather). The Mode-11 fold-spy (`test_L2_…`) still firing
    `2·n_levels` confirms the SAME fold path.
  - **Mode-12 lens:** `array_equal` on the bridged value pins the VALUE but is blind to CONTAINER TYPE —
    (ii)'s type assertions are the required companion (a FullField with equal flat values would pass (i)).
  - **Teeth:** (a) an interior↔boundary block swap in the re-type → `to_unified` reds (cells↔corner move);
    (b) a re-type minting Flux members instead of SourceSink → (ii) reds.

- **G-b3.2 — A_BA block applyᵀ ≡ old bulk pullback (P3).** `out = A_BA.apply_transpose(composite_cotangent)`
  (sphere, ≥2G, NONZERO seed cotangent — a present-zero seed gives `Reconstructionᵀ(0)=0`, blind):
  `array_equal(out.interior.values, _pullback_reconstruction(sn, S, chi_seed))` AND
  `out.radial_characteristic is None` (the forward-looking 2-block System-A shape).
  - **bit-id ruling:** BIT-IDENTICAL (the pullback `w·K_isoᵀ(Reconstructionᵀ χ)` is unchanged).
  - **Mode-12 lens:** `array_equal` pins the bulk pullback; the `radial_characteristic is None`
    assertion pins the new output shape (a transpose returning the old 3-slot shape would pass the
    value but leak a present-zero seed).
  - **Teeth:** (a) drop the pullback (present-zero bulk) → value reds; (b) return `radial_characteristic
    ≠ None` → the None-shape assertion reds.

- **G-b3.3 — adapter ≡ old operator byte-identical + THE SENTINEL-TEETH CHECK.** The FullField→FullField
  adapter (`_lagged_gains`' A_BA wrapper): `adapter.apply(psi)` (FullField→FullField) is `array_equal`
  to the OLD FullField output — `.radial_characteristic.values == _ba_oldloop_reference(emission, sn)`,
  `.interior`/`.boundary` present-zero. **THE SENTINEL-TEETH CHECK (confirms the L4-S reasoning):**
  wrap `RadialCharacteristicEmission.apply` (the CLASS method) with a counting spy, call `adapter.apply(psi)`
  ONCE in isolation, assert the counter == 1 — this PROVES the adapter DELEGATES to the wrapped instance's
  class method (not an inlined fold copy). This is the load-bearing assumption that makes the L4-S sentinel's
  spy fire through the adapter (Mode-11: a gate whose spy never reaches the rewired reader is vacuous).
  - **bit-id ruling:** BIT-IDENTICAL (`array_equal` on the whole FullField flat).
  - **Mode-12 lens:** `array_equal` on the full flat pins value AND placement; the spy-count is a Mode-11
    routing proof (blind to value — value is G-b3.1's job).
  - **Teeth:** (a) an adapter that DROPS the embed (present-zero ray) → `array_equal` reds; (b) an adapter
    that INLINES the fold instead of delegating → the spy-count check reds (counter stays 0).

- **G-b3.4 — B_b.H reciprocity through the NEW ray-composite space (Mode-12 CLOSED, RULING P2).**
  Re-express `test_euclidean_transpose_is_the_vcell_hilbert_adjoint` on System B's member space:
  `⟨B_b x, y⟩_G = ⟨x, B_bᵀ y⟩_G` under `G = sn.radial_characteristic_composite_space` (interior
  `V_cell` ⊕ boundary `V(r=R)`), `x, y` random RadialCharacteristicComposites. CONTROL `< 1e-12` +
  TWO teeth (a wrong-direction transpose; an asymmetric corner gauge) — both `> 1e-3`.
  - **bit-id ruling:** the reciprocity DEFECT is `< 1e-12` (Euclidean = Hilbert under the symmetric
    gauge, P2) — an equality claim, not a tolerance-relaxation.
  - **Mode-12 lens:** reciprocity is a SCALAR functional — Mode-12 CLOSED by the SPD composite metric
    (A2b) ONLY WITH the two teeth (which prove the functional is OUTSIDE the error class's stabiliser
    for the swap-direction and gauge-asymmetry mutations). Control-alone is fooled by a broken baseline;
    mutated-alone is fooled by a broken control — BOTH legs mandatory (L18).
  - **Teeth:** inherited verbatim from the existing gate (wrong-transpose = `_g_recip(fwd, fwd, g)`;
    asymmetric gauge = double the `+1` corner weight).

### Re-points (b3) — exactly how

- **`TestB_b_RayBoundary` (:485-583, all 6 rows).** B_b.apply/apply_transpose now take/return
  `RadialCharacteristicComposite`. Add a `_ray_composite(sn, seed_vals)` helper (System B with the
  seed placed in `interior`/`boundary`) replacing `_seed_composite`; a `_dense_ray_composite(fn, sn)`
  replacing `_dense_seed` (probe the composite flat layout). The corner-swap gate reads
  `out.boundary.corner(level, sign)`; the cells-zero gate reads `out.interior.cells(level, sign)`.
  The G_sd-reciprocity gate is subsumed by / re-expressed as G-b3.4 over the composite space.
- **`TestCoupledLift` L1-FWD / L1-ADJ / L2 / L3 (:1585-1833).** `_a_ba_scatter(sn).apply(psi)` now
  returns a `RadialCharacteristicComposite` — read `.interior` (the ray source) via `to_unified()`
  for the value oracle (`_ba_oldloop_reference`); the "A_BA bulk present-zero (no double-count)" check
  becomes STRUCTURAL (unspellable — the codomain has no bulk). L1-ADJ's `chi_cot` re-points from a
  FullField seed-only cotangent to a System-B composite cotangent; `A_BA.apply_transpose(chi_cot).interior`
  (bulk pullback) stays a FullField with `radial_characteristic=None`. The fwd↔adj Euclidean reciprocity
  (leg c) uses `A_BA.apply(psi).to_unified().values @ chi_seed`.
- **`TestCoupledLift` isinstance gates (:1848 sphere, :1858 slab) + the L4-S sentinel (:1837-1900).**
  After b3 the gain is the ADAPTER, so `any(isinstance(g, RadialCharacteristicEmission) …)` is now
  FALSE — re-point to a wraps-A_BA predicate (`isinstance(getattr(g, "_emission", g),
  RadialCharacteristicEmission)`, or `isinstance(g, <AdapterType>)`). The Mode-11 spy on
  `RadialCharacteristicEmission.apply` STILL FIRES through the adapter (class-method dispatch —
  see G-b3.3's proof) — the `counter > 0` leg is UNCHANGED and correct. The slab (:1858) leg stays a
  seedless control (no A_BA, no adapter). The L4-S teeth test (:1877, un-widened `_lagged_gains`)
  is UNCHANGED (drops the whole gain slot → counter 0).
- **`test_g_adjoint_reciprocity.py::test_tooth_a_ba_transpose_drop_reds (:660) + `_full_loss_case`.**
  The composite loss `A = (L+C) − S − A_BA − B` must remain composable over FullFieldSpace → the
  A_BA and B_b terms in the loss build are the ADAPTER-wrapped operators (FullField→FullField). The
  `_drop_pullback` monkeypatch on the bare `RadialCharacteristicEmission.apply_transpose` still fires
  through the adapter (same dispatch as L4-S); re-point its RETURN to the b3 shape (FullField,
  `radial_characteristic=None`). The A2a/A2b composite reciprocity gates likewise route the loss through
  the adapters (no change to the metric arithmetic — the adapters are behavior-neutral by construction).
- **`TestRegressionFloor::test_bulk_to_ray_coupling_lives_in_the_lagged_A_BA_gain (:214).** `_dense(A_BA.apply,
  tpl)` BREAKS — a FullField-templated probe cannot densify a FullField→RadialCharacteristicComposite map
  (`.to_flat()` size mismatch). Re-point to `_dense(adapter.apply, tpl)` (the FullField→FullField face):
  `A_BA_sb` (bulk→ray) is measured through the adapter; `A_BA_bb` (bulk→bulk = 0) is now STRUCTURAL (the
  adapter writes only the ray slot). `test_loss_operator_is_block_triangular_in_the_ray (:194)` is
  UNAFFECTED (it probes `(L+C)`, which stays FullField→FullField — the fused walk keeps the ray as System
  A's third block internally through 4e).

### UNAFFECTED (confirm, do not touch)

- **`TestA_BA_SchurFold` (:1307-1508).** Pins the fold FACTOR (`_apply_A_BA` =
  `RadialCharacteristicReconstruction.apply`, the `# BIND:` surface), NOT the full A_BA operator —
  unaffected by the operator re-type. CONFIRM `_apply_A_BA`/`_apply_A_BA_transpose` still route through
  the Reconstruction factor (they do — they construct `RadialCharacteristicReconstruction(sn)` directly).
- **B.1d `TestSplitFidelityBridge` Flux round-trips, B.1c flux/displacement leaf laws.** Additive; stay
  green byte-identical.

---

## Mutation ledger (monkeypatch-only; NEVER `git checkout` an uncommitted file — L4)

| sub | gate | tooth (in-process) | target `file:line` (⚠ best-estimate) | expected RED |
|---|---|---|---|---|
| b1 | G-b1.1 | SourceSink inherits FluxRole | NEW source-sink leaf class | source+source positive control reds |
| b1 | G-b1.3 (CRUX) | role-table SourceSink→Flux | `radial_characteristic_composite.py:122` (`from_unified`) | SourceSink/Displacement ROLE reds; Flux GREEN (asymmetry = evidence) |
| b1 | G-b1.3 | bridge drops the corner | `from_unified` corner loop | `to_unified` boundary round-trip reds |
| b1 | G-b1.4 | revert static param to Flux leaf | `radial_characteristic_composite.py:70-72` | pyright `reportArgumentType` on source-member construction |
| b2 | G-b2.1 | presence-dispatch inversion (always-pass kwarg) | `full_field_space.py:328/336` (`_rebuild`) | `apply_metric(2-block)` raises |
| b2 | G-b2.1 | drop seed kwarg for FullField too | `_rebuild` presence branch | FullField ψ½-block metric reds |
| b2 | G-b2.2 | drop boundary block from shape sum | `full_field_space.py:239-248` (`from_blocks`) | coextensiveness reds |
| b3 | G-b3.1 | interior↔boundary block swap in re-type | `radial_characteristic.py:1276` (`RadialCharacteristicEmission.apply`) | `to_unified` value reds |
| b3 | G-b3.2 | drop the bulk pullback | `radial_characteristic.py:1335` (`apply_transpose`) | pullback value reds |
| b3 | G-b3.3 | adapter drops the embed / inlines the fold | NEW adapter in `solver.py:589-620` | `array_equal` reds / spy-count 0 |
| b3 | G-b3.4 | wrong-transpose + asymmetric gauge | inherited (`_g_recip`) | recip defect `> 1e-3` |

## Refutations (fire before the ink dries)

1. **B.2b is a RE-LABELING — `array_equal`, never `rtol`/`nulp`.** No arithmetic moves; a tolerance
   in a B.2b value gate would hide a re-labeling bug that DID move arithmetic. Principled-equiv enters
   at B.2d/4e, not here.
2. **G-b1.3 role-preservation is (role class) ⊕ (values), each blind to the other.** The current hardcoded
   bridge passes the value round-trip for ALL roles yet is WRONG for SourceSink/Displacement — the ROLE
   assertion is the ONLY catcher, and its RED-on-SourceSink / GREEN-on-Flux asymmetry is the config-blindness
   evidence (§0.6). A Flux-only bridge test is structurally blind.
3. **The b2 `_rebuild` presence-dispatch is a BRANCH — exercise BOTH arms.** A suite testing only the FullField
   (seed) arm OR only the 2-block arm is blind to the wrong branch (L11); the metric-trio gate MUST run on a
   REAL RadialCharacteristicComposite (2-block) AND assert FullField byte-identity (seed).
4. **The L4-S spy fires through the adapter ONLY IF the adapter delegates to the wrapped class method.**
   G-b3.3's isolated spy-count == 1 PROVES delegation (not inline); without it the L4-S sentinel is
   Mode-11-vacuous for the NEW adapter path.
5. **`_dense(A_BA.apply, tpl)` breaks after the codomain change** — re-point the RegressionFloor A_BA probe
   to `_dense(adapter.apply, tpl)`; `A_BA_bb=0` becomes STRUCTURAL, not a measured near-zero.
6. **cyl/slab are the non-carrying CONTROL** — the SourceSink leaves / composite / member-space are
   UNCONSTRUCTABLE there (the presence-raise IS the gate). All value rows ride the carrying sphere-GL S4.
7. **The re-type makes "A_BA writes bulk" / "B_b touches trace" UNSPELLABLE (Pattern-4)** — the old
   present-zero-disjointness gates become structural type facts; assert the container type, do not measure
   a near-zero bulk block that no longer has a slot.
8. **`_ba_oldloop_reference` is a bit-id INHERITANCE oracle (shares the fold kernel), sufficient for a
   re-labeling** — it catches the container/placement bug (the only thing the re-type can break). The
   end-to-end structural anchor (φ=Q/Σ_t, k_inf) is the existing suite's wall, not re-proven per gate.

## Result contract

Land in three sub-commits (b1 leaves+bridge / b2 member space / b3 operator re-type + adapters). NEW test
homes: extend `tests/sn/mesh/test_radial_characteristic_split_leaves.py` (G-b1.1/2), `tests/transport/
test_radial_characteristic_composite.py` (G-b1.3/4), NEW `tests/sn/mesh/test_radial_characteristic_composite_space.py`
(G-b2.1/2); extend `tests/sn/operators/test_psi_half_coupling.py` (G-b3.1-4 + all the re-points) +
`tests/sn/operators/test_g_adjoint_reciprocity.py` (`test_tooth_a_ba_transpose_drop_reds` re-point). Every
tooth mutation-verified in-process under `-O`. End-to-end acceptance: full `tests/sn -m "not slow"` +
`tests/numerics` GREEN + ratchet `transport:1` + pyright 0 on the new source-sink leaves (G-b1.4 static) +
sphinx -W. The load-bearing deliverables: **G-b1.3 role-preservation (role⊕values split)** + **G-b2.1
presence-dispatch metric-trio (both branches)** + **G-b3.3 the adapter-delegation sentinel-teeth check**
(the proof that keeps the L4-S sentinel non-vacuous).
