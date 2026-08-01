# B3.2 — landing the narrowed contract (test migration)

Execution memo for campaign phase **B3.2** (SN boundary law narrowed to
`Γ₊ → Γ₋`). Branch `refactor/operator-strategy-layers`, production at
`b39502f8` + uncommitted `orpheus/` edits. Repo root
`/Users/rodrigo/git/nuclear/ORPHEUS`, `.venv/bin/python -O -m pytest`, SERIAL.

Design of record: `.claude/plans/b3_domain_narrowing_crosswalk.md`.
Gate design: `scratch/b3_gate_reposing.md` (RG-1…RG-5). Lessons: L29.

Evidence key: **[M]** measured this session (pytest / probe output pasted) ·
**[R]** read (`file:line`) · **[G]** grep.

**Scope discipline.** Nothing under `orpheus/` written. No `git checkout` /
`restore` / `stash` / `clean` on any path. Mutations are in-process
monkeypatches or tmp copies.

---

## (c) GENUINE VALUE REGRESSIONS — **NONE FOUND**

No failure in the inherited 49 is a value regression. Every one is a
signature/API break or a stale expectation. Positive evidence rather than
absence-of-evidence:

* **[M]** The wide baseline over `tests/{geometry,sn,diffusion,numerics,transport}
  -m "not slow"` reproduces **exactly the same 49** —
  `49 failed, 4251 passed, 6 skipped, 111 deselected, 57 xfailed … in 729.02s`.
  Nothing beyond the boundary-construction tests moved: **zero** solver,
  adjoint-reciprocity, DSA, Gauss-Seidel, snapshot-regression or convergence
  failures.
* **[M]** All 49 reds are confined to *realizer-construction* and
  *law-object-probing* tests.
* **[M]** The composite `SNBoundaryOperator.apply` / `.apply_transpose` values
  are byte-identical — proven in §3 against a reference materialised from the
  **retired** expression, not by calling the new code twice.

**[M] FINAL STATE** (same scope, after the migration):

```
4335 passed, 6 skipped, 111 deselected, 73 xfailed, 21 warnings in 985.85s (0:16:25)
```

`4251 → 4335` passed = 49 re-posed + 35 new. `57 → 73` xfailed = +16, exactly
the deliberate strict-xfail set (4 RG-3b + 7 domain gate + 5 mixed-BC).

Three **findings** that are *not* (c) but that the phase owes an answer to are
in §5. The load-bearing one: **the narrowed law does not validate its domain**
— a full-face input is silently accepted and returns wrong values.

---

## 0. The inherited state, measured

**[M]** `python -O -m pytest tests/sn/operators tests/geometry -q -p no:randomly -m "not slow"`:

```
49 failed, 1130 passed, 4 skipped, 1 deselected, 2 xfailed, 27 warnings in 18.63s
```

Reproduces the reported number exactly.

---

## 1. Classification — every failure, one bucket

Classification rule (decidable, applied uniformly): **(a)** iff the test passes
*unchanged* once the construction call is migrated; **(b)** iff the *assertion*
must change to state the narrowed contract. (b) dominates when both symptoms
are present.

### 1.1 Bucket (a) — signature/API migration, 4 items

| # | test | why | fix |
|---|---|---|---|
| a1–a3 | `test_operator_block_role.py::TestBoundaryLeaves::test_realized_bc_advertises_boundary_operator[vacuum, reflective_albedo1, reflective_albedo07]` | `_boundary_method_space()` builds `SNMethodSpace(quadrature, face, inflow_indices)` — no `outflow_indices` | add `outflow_indices` to the helper. Assertion (`block_role is BOUNDARY`) unchanged |
| a4 | `test_bc_errors.py::test_boundary_error_base_class_contract` | its **positive control** realizes `ReflectiveBoundary` on `SNMethodSpace.minimal(quad)` | migrate the control to a face-ful space; the three `raises` legs unchanged |

### 1.2 Bucket (b) — stale expectation, re-posed, 45 items

Grouped by the retired representation each asserts.

**(b-I) "vacuum preserves the outflow rows"** — the `IncomingOrdinateMaskTensor`
pass-through contract. Retired: vacuum is now the zero map `Γ₊ → Γ₋`.

| test | old assertion | re-posed as |
|---|---|---|
| `test_sn_boundary_realizer.py::TestRealizeVacuum::test_lebedev17_xmin_zeroes_only_inflow` | `out[non_inflow] == psi[non_inflow]` | image is `zeros((|Γ₋|,)+tail)`; **and** `apply(full)` must not emit `N` rows |
| `…::test_lebedev17_ymax_zeroes_only_inflow` | same | same |
| `…::test_vacuum_returns_tensor_product` | `isinstance(ops[0], IncomingOrdinateMaskTensor)` | `isinstance(op, ZeroOperator)` + the two space hooks emit `Γ₋` / `Γ₊` |
| `test_bc_equivalence_snapshot.py::TestVacuumLebedev17Snapshot::test_realizer_zeroes_only_inflow_per_section_16A5` | outflow/tangential pass through | narrowed image `== 0`, shape `(|Γ₋|,…)`; the snapshot's own `psi_out` still supplies the input |
| `test_boundary.py::test_vacuum_bc_realizer_zeros_only_inflow_per_section_16A5` | `psi_in[outflow] == psi_out[outflow]` | same |
| `test_snmesh_realizer_wiring.py::test_2d_cartesian_vacuum_xmin_masks_only_inflow` | `out[non_inflow] == psi[non_inflow]` | same, through the mesh-wired shim |
| `…::test_1d_cartesian_vacuum_right_masks_only_inflow` | same | same |
| `…::test_1d_spherical_vacuum_routes_through_realizer` | `_angular_factor(inner)` is `IncomingOrdinateMaskTensor` + pass-through | `inner` is `ZeroOperator`; narrowed image is zero |
| `test_bound_compat.py::Test188WiringContracts::test_curvilinear_realize_boundary_law_routes_through_realizer` | `inner.ops[0]` is `IncomingOrdinateMaskTensor` | `inner` is `ZeroOperator` |
| `test_boundary.py::test_albedo_zero_and_vacuum_agree_on_inflow_rows` | both zero the inflow rows of a full face | both zero on `Γ₋`, each on its own domain |

**(b-II) "reflective is a full-`N` permutation"** — retired: the perm is on the
REDUCED ordinate axis (`local_perm = γ₊.to_local(perm[inflow])`).

| test | old assertion | re-posed as |
|---|---|---|
| `test_sn_boundary_realizer.py::TestRealizeReflective::test_specular_unit_albedo_lebedev_matches_hand_computed` | `op.apply(psi) == psi[reflection_index]` | `op.apply(γ₊ψ) == psi[perm][inflow]` — the **retired expression, restricted**; independent of `to_local` |
| `…::test_specular_partial_albedo_matches_hand_computed` | `0.7 * psi[perm]` | `0.7 * psi[perm][inflow]` |
| `…::test_specular_unit_albedo_returns_tensor_product` | TP(`PermutationOperator`, `Identity`) | unchanged type; **plus** `perm.size == |Γ₊|` (the narrowing is now type-visible) |
| `…::test_specular_partial_albedo_returns_scaled_tensor_product` | `ScaledOperator(0.5, TP)` | same + reduced-perm size |
| `test_snmesh_realizer_wiring.py::test_2d_cartesian_reflective_ymax_returns_permutation` | `apply(psi) == psi[ref]` | `apply(γ₊ψ) == psi[ref][inflow]` |
| `…::test_2d_reflective_y_face_builds_y_axis_permutation` | `perm.perm == reflection_index(axis)` | `perm.perm == γ₊.to_local(reflection_index(axis)[inflow])`, with the x≠y non-vacuity guard **kept and strengthened** |
| `…::test_1d_cylindrical_one_boundary_outer_reflective` | `isinstance(angular, PermutationOperator)` | + reduced size |
| `test_bound_compat.py::…::test_1d_reflective_faces_realized_as_permutation_tp` | `angular.perm == reflection_index("x")` | reduced perm; **slab is the discriminating fixture for `to_local`** (see §3) |
| `test_boundary.py::test_specular_bc_indexes_through_reflection_partner` | `psi_in == psi_out[ref]` | narrowed |
| `…::test_specular_bc_with_partial_albedo` | `0.5 * psi_out[::-1]` | narrowed |
| `…::test_specular_bc_axis_y_on_lebedev` | `psi_out[ref_y]` | narrowed |
| `…::test_specular_realized_op_advertises_apply_transpose` | `realized.is_adjointable` | unchanged assertion, face-ful space |
| `…::test_specular_apply_transpose_reciprocity_unweighted` | `⟨Bψ,φ⟩ == ⟨ψ,Bᵀφ⟩` on a full face | the **rectangular** reciprocity: `⟨B γ₊ψ, φ_Γ₋⟩ == ⟨γ₊ψ, Bᵀφ_Γ₋⟩` |
| `…::test_specular_self_inverse_identity` | `B(B(x)) == α²x` | **RETIRED CLAIM, re-posed**: a `Γ₊→Γ₋` map is not composable with itself. The involution now lives on `γ₊ᵀ∘B` vs the mirror — see §2.3 |
| `test_bc_equivalence_snapshot.py::TestSpecularXLebedev17Snapshot::test_realizer_apply_matches_snapshot` | full-face `psi_in` | `snapshot["psi_in"][inflow]` — the frozen pre-B3.2 baseline, restricted |
| `…::TestSpecularYPartial07LS6Snapshot::test_realizer_apply_matches_snapshot` | same | same |
| `test_law_composition.py::test_realize_recursively_leaf_dispatches_to_sn_realizer` | walker ≡ direct on a full face | on `Γ₊` |
| `…::test_realize_recursively_lawscaled_wraps_in_scaled_operator` | ditto | on `Γ₊` |

**(b-III) blocked on B3.4 — a narrowed law cannot compose with an un-narrowed
one.** 6 items. **[M]** `0.3·spec + 0.7·white` now raises on BOTH domains:

```
0.3*spec+0.7*white on Gamma+(4, 2) -> RAISE ValueError: AngularAverageOperator.apply: psi.shape[0] = 4, expected 8
0.3*spec+0.7*white on FULL(8, 2)   -> RAISE ValueError: operands could not be broadcast together with shapes (4,2) (8,2)
```

| test | disposition |
|---|---|
| `test_law_composition.py::test_realize_recursively_lawsum_returns_operator_sum` | structural-only (no `apply`) → **migrates**, stays green |
| `…::test_realize_recursively_walks_nested_depth_first` | structural-only → **migrates** |
| `…::test_realize_recursively_apply_matches_pointwise_weighted_sum` | `xfail(strict=True)` naming **B3.4 / #183 / #189** |
| `…::test_realize_recursively_nested_apply_matches_distributive_form` | `xfail(strict=True)` naming B3.4 |
| `test_boundary.py::test_wave0_sum_of_realized_bcs_acts_as_weighted_sum` | `xfail(strict=True)` naming B3.4 |
| `…::test_operator_sum_of_bcs_acts_as_weighted_sum` | `xfail(strict=True)` naming B3.4 |
| `test_bc_equivalence_snapshot.py::TestMixed30Spec70WhiteLS4Snapshot::test_realizer_apply_matches_snapshot` | `xfail(strict=True)` naming B3.4 |

**(b-IV) the vacuum orientation guard's two positive controls** —
`test_sn_boundary_realizer.py::TestVacuumTraceOrientationGuard::{test_correct_orientation_realizes,
test_faceless_space_carries_no_orientation_truth}`. Both assert
`isinstance(op, TensorProductOperator)`; retired → `ZeroOperator`. The
*faceless-escape* claim survives verbatim (a hand-built space with BOTH index
sets and no `face` still realizes, guard silent) — see §2.4.

**(b-V) `test_apply_per_face_equals_legacy_bc_apply_on_inflow_row` ×4 and
`test_adjointable_when_all_faces_support` ×4** — the C-1 gates proper. Re-posed
as RG-1/RG-2/RG-2b and G3′/RG-4. §2.

---

### 1.3 Where each item landed

| file | before | after |
|---|---|---|
| `tests/sn/operators/test_sn_boundary_operator.py` | 8 red | RG-1/RG-2/RG-2b/G3′/RG-4/RG-3/RG-3b/RG-5 — **38 passed, 4 xfailed** |
| `tests/sn/operators/test_sn_boundary_realizer.py` | 9 red | **26 passed** (+2 NEW negatives) |
| `tests/sn/operators/test_snmesh_realizer_wiring.py` | 6 red | **11 passed** |
| `tests/sn/operators/test_operator_block_role.py` | 3 red | **20 passed** |
| `tests/geometry/test_boundary.py` | 10 red | **23 passed, 2 xfailed** (B3.4) |
| `tests/geometry/test_law_composition.py` | 6 red | **16 passed, 2 xfailed** (B3.4) |
| `tests/geometry/test_bc_equivalence_snapshot.py` | 4 red | **7 passed, 1 xfailed** (B3.4) |
| `tests/geometry/test_bound_compat.py` | 2 red | **13 passed** |
| `tests/geometry/test_bc_errors.py` | 1 red | **11 passed** |
| `tests/sn/operators/test_b3_domain_narrowing.py` | — | **NEW: 20 passed, 7 xfailed** |

**Nothing deleted. Nothing skipped. No tolerance touched** — every re-posed
assertion is `array_equal` or an exact predicate, exactly as before.

New shared helpers in `tests/sn/_test_helpers.py` (one source for the migrated
call sites — Pattern 2): `face_trace`, `face_method_space`, `local_positions`.

---

## 2. The re-posed contracts, and the mutation that reddens each

Harness: `scratch/b3_2_mutations.py` — an in-process pytest plugin, IN THE
REPO (the B3 plan's "the sweep still catches 30/31" gate was previously
unrunnable because the harness lived in job scratch and evaporated). Every
variant prints a sha256 bite-check of the composite's forward + transpose
output on a fixed fixture, so a patch that silently fails to install is
visible rather than being misread as "the gate is blind" (the vv Mode-8
METHOD WARNING).

```
ORPHEUS_B32=N1 .venv/bin/python -O -m pytest tests/sn/operators \
    -p no:randomly -p scratch.b3_2_mutations -q
```

### 2.1 Measured colour matrix — `test_sn_boundary_operator.py`

**[M]** CONTROL (no mutation): `38 passed, 4 xfailed, 1 warning in 0.17s`.

| mutation | items red | which gates |
|---|---|---|
| **N1** write the image to the OUTFLOW rows | 14 | RG-1 ×4, RG-2 ×4, RG-2b ×4, `test_subset_reflects_only_selected_faces`, RG-5 |
| **N2** pass `γ₊` through onto the outflow rows | 5 | RG-2 ×4, RG-5 |
| **N3** hand the law the WRONG half (`γ₋` for `γ₊`) | 4 | RG-1 ×4 |
| **N4** hand the law the FULL face | 4 | RG-1 ×4 |
| **N5** leak onto the TANGENTIAL rows | **1** | **RG-2[cyl] ONLY** |
| **N7** transpose writes onto `Γ₋` | 8 | G3′ ×4, RG-4 ×4 |
| **N8** split remap → `arange(sel.size)` | **1** | **RG-5 ONLY** |
| **M1/M2** reflective permutation content | **0** | — see §2.4 |

### 2.2 ⭐ RG-2's widening is MEASURABLY necessary, not cosmetic

The pre-B3.2 leg was `assert not got[outflow].any()`. It survives B3.2 and can
still red — but it is blind to the tangential rows. **[M]** on the real
production shape (`/tmp/b32/probe_n5.py`, `cyl_reflective`, `product(2,4)`):

```
  control  face=xmax  |I|=2 |O|=2 |T|=4
      OLD leg  got[outflow].any() -> False   (False == GREEN)
      NEW leg  max|got[complement(inflow)]| -> 0.000e+00   (0 == GREEN)
  N5       face=xmax  |I|=2 |O|=2 |T|=4
      OLD leg  got[outflow].any() -> False   (False == GREEN)
      NEW leg  max|got[complement(inflow)]| -> 1.846e+00   (0 == GREEN)
```

A 1.846 leak onto `Γ`'s tangential band is **GREEN under the old spelling**.
`outflow` → `complement(inflow)` costs nothing and closes it. The face slot is
`I ⊔ O ⊔ T`, and "not inflow" must never be spelled "outflow".

### 2.3 ⭐ Both `to_local` sites are covered, and neither fixture covers the other

The naive `arange(sel.size)` appears at TWO sites. Measured across
`tests/{sn/operators,geometry} -m "not slow"`:

| mutation | items red | verdict |
|---|---|---|
| **N8** — SPLIT remap → `arange` | **1** | **only** `TestScheduleSplitPartition2D` (RG-5). The pre-existing sphere-based split gate `test_psi_half_coupling.py::test_split_masked_halves_are_trace_only` stays GREEN — 1-D is blind, as the crosswalk predicted |
| **N9** — REFLECTIVE remap → `arange` | **16** | the SLAB rows red (`test_bound_compat.py::test_1d_reflective_faces_realized_as_permutation_tp`, `test_boundary.py::test_specular_bc_with_partial_albedo`, both slab bit-identity rows); **`cyl_reflective` and `cart2d_mixed` bit-identity stay GREEN** — those quadratures are exactly where `perm[inflow]` is already `outflow` in order |

So: the **slab** discriminates the reflective remap and is blind to the split
remap; **2-D** discriminates the split remap and is blind to the reflective
one. Both gates ship, and each carries an ACTIVATION guard (`discriminating`
counter in `test_bound_compat.py`; `test_the_2d_lower_rows_are_not_a_prefix_of_the_inflow_set`
beside RG-5) so neither can silently decay into the other's blind fixture.

### 2.4 The honest division of labour — M1/M2

**M1** (reflective perm → identity) and **M2** (rolled by one) BIT
(`FWD sha 9f3e6044b46abe3b → 3563a9c81be3de48`) but redden **nothing** in
`test_sn_boundary_operator.py`. That is by design, not a hole: RG-1's
reference is `bc.apply(γ₊ψ)`, so it gates the WIRING (which domain, which
rows) and takes the law's math as given. The law's math is gated by the
hand-computed references in the realizer/wiring/snapshot files and by the
bit-identity module — **[M]** M1 reds 6 and M2 reds 10 items in
`test_b3_domain_narrowing.py`. Stated here so a future reader does not read
RG-1's green under M1 as coverage.

### 2.5 Re-posed contracts, one line each

| gate | the narrowed contract | reddens on |
|---|---|---|
| **RG-1** `test_apply_per_face_equals_law_of_the_outflow_half_trace` | `B[inflow] == law(γ₊ψ)` | N1, N3, N4, face↔face swap |
| **RG-2** `test_apply_emits_nothing_outside_the_inflow_trace` | `B[complement(inflow)] == 0` via `pytest.fail` | N1, N2, **N5** |
| **RG-2b** `test_apply_actually_emits_on_the_inflow_trace` | `Γ₋` is not all-zero | N1 |
| **RG-3** `test_realized_law_maps_gamma_plus_to_gamma_minus` | `(|Γ₊|,…) → (|Γ₋|,…)` **and** not an endomorphism | a revert to the full-`N` law |
| **RG-3b** `test_realized_law_refuses_a_full_face_input` | full-face input RAISES | **strict xfail — the guard is MISSING, §5.1** |
| **G3′** `test_adjointable_when_all_faces_support` | `Bᵀ == ι₊ ∘ lawᵀ ∘ γ₋` | N7 |
| **RG-4** `test_transpose_emits_nothing_outside_the_outflow_trace` | `Bᵀ[complement(outflow)] == 0` | N7 |
| **RG-5** `test_split_halves_partition_the_whole_trace_on_a_2d_mesh` | `B_lower + B_upper == B`, 2-D | **N8 (sole catcher)** |
| **bit-id** `test_composite_is_bit_identical_to_pre_b32` | `array_equal` vs the RETIRED expression | N1–N5, N7, N9, **M1, M2** |

### 2.6 The `_NoTransposeLaw` duck-typed surrogate (G4)

`tests/sn/operators/test_sn_boundary_operator.py:…::_NoTransposeLaw.apply` was
a full-face identity (`return x`). It never failed — the transpose path raises
before the body runs — so it is invisible to every search the retirement audit
performs (grep, call graph). Narrowed to emit `zeros((|Γ₋|,) + tail)`, with the
reason inline. This is the "fourth search" class the B2 audit added.

---

## 3. Bit-identity evidence

`tests/sn/operators/test_b3_domain_narrowing.py::TestBitIdentityAgainstTheRetiredExpression`.

* **Reference**: `_pre_b32_face_action` — the retired expression transcribed in
  numpy against the law DESCRIPTOR (`sn.bc[face].law`). Vacuum = "zero the
  inflow rows, keep the rest, then write only the inflow rows" ⟹ identically
  zero; reflective = `α · np.take(face_in, reflection_index(axis), 0)` then
  `[inflow]`. The transpose reference uses `np.argsort(perm)`, **not**
  production's cached `inverse_perm`. Nothing B3.2 touched is on the reference
  side.
* **Contract**: `np.array_equal`. The narrowing removes ROWS from a gather and
  changes no reduction tree, so the predicted drift bound is exactly zero —
  a `nulp` tolerance here would hide the class the gate exists for.
* **Fixtures**: slab (vacuum+reflective, asymmetric), slab (reflective ×2),
  sphere, **cylinder** (`product(2,4)` — 4 tangential of 8), **2-D Cartesian**
  (`level_symmetric(4)`, `ng=2`, mixed vacuum/reflective, `nx=4 ≠ ny=3`).
  Both directions. 10 rows.
* **Result [M]**: `20 passed, 7 xfailed, 1 warning in 0.15s` — every row green.
* **Non-vacuity**: each row asserts the reference is not identically zero
  across the faces, and a companion `test_the_retired_reference_is_falsifiable`
  perturbs `Γ₊` and requires the composite to move.
* **Teeth [M]**: reds under N1 (5), N2 (5), N3 (10), N4 (7), N5 (1), N7 (5),
  N9 (6), M1 (6), M2 (10).
* **Independent anchors, still green**: the closed-form gates downstream
  (`keff == k_inf`, `φ == Q/Σ_t`) never moved — old-vs-new is necessary, not
  sufficient (lesson L2), and these are the structurally-independent side.

---

## 4. The seven-law domain gate — B3.4's todo list

`tests/sn/operators/test_b3_domain_narrowing.py::TestEveryLawsDomain`.

Two legs per law. **Leg A**: `apply` on a `|Γ₊|`-row probe returns `|Γ₋|` rows.
**Leg B (the anti-Mode-12 leg)**: fed a FULL face, the law must not emit `N`
rows. Leg B is load-bearing because **[M]** `|Γ₊| == |Γ₋|` on EVERY quadrature
× face pair in the tree — gauss_legendre 4/5/8, product 2×4/3×4/4×8, lebedev
9/17, level_symmetric 4/6, all faces — so a shape assertion alone cannot tell
`Γ₊ → Γ₋` from `Γ₊ → Γ₊`. Leg B is exactly what catches `albedo` and
`periodic`, which pass Leg A by accident.

### 4.1 ⚠ REFUTATION of the brief's premise: it is SIX rows across FOUR law kinds, not three

The brief names "white, albedo(α∉{0,1}), periodic". Measured:

| law | narrowed? | measured behaviour on `Γ₊` / on the FULL face |
|---|---|---|
| `vacuum` | ✅ | `ZeroOperator` with both space hooks → `(|Γ₋|,…)` |
| `reflective` α=1, α=0.7 | ✅ | reduced `PermutationOperator` → `(|Γ₋|,…)` |
| `white` α=1, α=0.3 | ❌ | **RAISES** `AngularAverageOperator.apply: psi.shape[0] = 4, expected 8` |
| `albedo` **α=0.0** | ❌ | bare `ZeroOperator` (no hooks) → echoes; FULL → `(8,…)` |
| `albedo` **α=1.0** | ❌ | `IdentityOperator` → echoes; FULL → `(8,…)` |
| `albedo` α=0.5 | ❌ | `ScaledOperator(0.5, I&I)` → echoes; FULL → `(8,…)` |
| `periodic` | ❌ | `PeriodicWrapOperator & I` → echoes; FULL → `(8,…)` |
| `prescribed_inflow` | ❌ | **RAISES** `IncomingSourceOperator: ordinate axis mismatch — the inflow mask covers 8 ordinates but the requested block has axis-0 length 4` |
| `zero_flux` | n/a | SN REFUSES it (`BoundaryError(law="zero_flux")`) — pinned as its own negative, not an xfail |

**α=0 and α=1 albedo are un-narrowed too** — the fast paths return a bare
`ZeroOperator` / `IdentityOperator`, both endomorphisms. And
**`prescribed_inflow` is a fourth un-narrowed kind** the brief did not name.
The xfail reason string states all of this.

### 4.2 The markers are proven to FLIP

**[M]** With a landing-simulation plugin that narrows ONLY albedo (returning a
reduced `PermutationOperator` / a two-hook `ZeroOperator` at α=0):

```
3 failed, 3 passed, 17 deselected, 4 xfailed, 1 warning in 0.12s
```

All three albedo rows become `XPASS(strict)`; white / periodic /
prescribed stay `XFAIL`. So the marker set is a live todo list that flips
exactly for the law B3.4 narrows and for no other.

A simulation artefact worth recording: the first attempt used
`ScaledOperator(0.0, base)` for albedo(0) and the row stayed xfailed — because
`ScaledOperator` **refuses a zero scalar** ("degenerate; use ZeroOperator
explicitly"). The xfail swallowed that, which is the Mode-8 class-4 trap
appearing inside the verification of a class-4 defence. Fixed by simulating
the honest two-hook `ZeroOperator`.

### 4.3 Completeness

`test_the_law_inventory_is_complete` cross-checks the gate's law list against
`BoundaryTraceLaw.registry`, so a law added tomorrow cannot escape the gate.

---

## 5. Findings the phase owes an answer to (NOT (c) — no value moved)

### 5.1 ⚠ The narrowed law does NOT validate its domain

**[M]** With `gauss_legendre(8)` at `xmax` (`|Γ₊| = |Γ₋| = 4`, `N = 8`):

```
vacuum      ZeroOperator            on-Gamma+(4, 2)->(4, 2)   on-FULL(N=8)->(4, 2)
reflective  TensorProductOperator   on-Gamma+(4, 2)->(4, 2)   on-FULL(N=8)->(4, 2)
```

Both **silently accept a full-face input and return `|Γ₋|` rows of WRONG
values, with no raise**. `TraceRestrictionOperator` carries the shape guard the
crosswalk §9 designed — but the operator the realizer EMITS is a bare
`PermutationOperator` on the reduced axis (`np.take` truncates) or a
`ZeroOperator` whose hook ignores its input length.

Not a (c): the composite always composes `γ₊` first, so nothing production
computes has moved. But the crosswalk's own headline mitigation ("B3's own
headline bug is silent unless the carve adds a guard") is **not yet in the
tree**. Shipped as **RG-3b**, `xfail(strict=True)` ×4, verified to red for its
documented reason (`DID NOT RAISE any of (ValueError, IndexError)`) at the
`pytest.raises` line. Delete the marker when the guard lands.

### 5.2 A narrowed law can no longer compose with an un-narrowed one

**[M]** `0.3·spec + 0.7·white` raises on BOTH domains. 6 mixed-BC gates are
therefore B3.4-blocked and ship as strict xfails naming it (each verified with
`--runxfail` to red at the documented `AngularAverageOperator.apply` refusal,
after re-posing their probes onto `Γ₊` so the failure is attributable rather
than a reference-side broadcast error). This is a real architectural
consequence of narrowing only two of seven laws — worth a ruling on whether
B3.4 should be sequenced immediately after B3.2 rather than after B3.3.

### 5.3 Stale prose in production (doc-only, no code change requested)

`orpheus/sn/operators/boundary.py::_reflect_trace`'s docstring still describes
the pre-B3.2 world: *"The realized per-face law … is a **full-face** operator,
so a non-zero outflow emission would corrupt …"*, and the whole
``B_face = P_inflow ∘ law`` / ``B_faceᵀ = lawᵀ ∘ P_inflow`` paragraph. The
inline comments in the body are correct and current; the docstring above them
is not. Flagged, not touched (no `orpheus/` edits in this pass).

`test_snmesh_realizer_wiring.py`'s module docstring likewise still says the
suite pins "Wave-2 trace-mask construction"; the mask is retired.

### 5.4 `SNMethodSpace.minimal` is now a partial constructor

Post-B3.2 it can realize white / albedo / periodic / prescribed but NOT vacuum
or reflective — the two laws SN actually reaches from a mesh. Its docstring
still says "Vacuum realization … raises BoundaryError … that's the intended
behaviour" and does not mention reflective. Once B3.4 narrows the rest,
`minimal` cannot realize ANY law and becomes a retirement candidate. Recorded
so B3.4 does not have to rediscover it.

---

## 6. Standing evidence

* Mutation harness: `scratch/b3_2_mutations.py` (in the repo — finding #7 of
  the gate-reposing memo discharged for this phase's gates).
* Probes (job scratch, reproducible from the harness): `/tmp/b32/probe1.py`
  (per-law domain behaviour), `/tmp/b32/probe2.py` (the `|Γ₊|` vs `|Γ₋|` hunt +
  composition), `/tmp/b32/probe_n5.py` (old-vs-widened outflow leg).


---

## 7. Final measured state

**[M]** `python -O -m pytest tests/{geometry,sn,diffusion,numerics,transport} -q -p no:randomly -m "not slow"`:

```
4335 passed, 6 skipped, 111 deselected, 73 xfailed, 21 warnings in 985.85s (0:16:25)
```

**[M]** `python -O -m pytest tests/sn/operators tests/geometry -q -p no:randomly -m "not slow"`:

```
1214 passed, 4 skipped, 1 deselected, 18 xfailed, 27 warnings in 23.61s
```

**[M]** `python -O -m tests._harness.audit` — `Orphan equations (0 of 325
testable theory labels have zero test coverage)`, `error_catalog.md ERR
coverage (71/71 entries have a catching test)`. No V&V regression; the new
module is `foundation`-tagged and carries no `verifies(...)` (lesson L9 — a
structural invariant with no theory-page `:label:` must not stack both).

### The 16 deliberate strict xfails — the todo list B3.4 / the domain guard own

| count | where | deleted when |
|---|---|---|
| 4 | `test_sn_boundary_operator.py::TestNarrowedLawDomain::test_realized_law_refuses_a_full_face_input` | the narrowed law gets a domain guard (§5.1) |
| 7 | `test_b3_domain_narrowing.py::TestEveryLawsDomain` — white ×2, albedo ×3, periodic, prescribed | **B3.4** narrows each law |
| 2 | `test_boundary.py` mixed-BC linearity | B3.4 |
| 2 | `test_law_composition.py` walker distributivity | B3.4 |
| 1 | `test_bc_equivalence_snapshot.py::TestMixed30Spec70WhiteLS4Snapshot` | B3.4 |

Every one verified with `--runxfail` to red for ITS documented reason, and the
B3.4 set additionally proven to FLIP under a landing simulation (§4.2).
