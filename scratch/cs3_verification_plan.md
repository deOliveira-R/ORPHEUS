# CS3 verification plan — flux lives in the cone, and every gate that says so can RED

**Written PRE-carve** (proactive `test-architect` dispatch, `delegation.md` operator-algebra
trigger) for phase **CS3** of `.claude/plans/space_and_kernel_binding_campaign.md` §4.
Branch `refactor/cone-field-algebra`, HEAD **`000cf144`**, 2026-08-19.

Evidence base: `scratch/omr_v2_grounding/D_flux_ontology.md` (the blast-radius memo). Every
`[M]` below is MINE, re-measured at `000cf144` — the memo's structure held; where I widened
or refuted it the row says so.

**How to read.** §1 is the value-neutrality harness (what already pins each production
consumer, and what does not). §2 is the ρ-trajectory capture gate — **already written,
committed to the tree, green** (`tests/numerics/test_si_diagnostic_trajectory.py`). §3 is
the per-file migration map. §4 is the new-algebra gate set. §5 is the cone predicate. §6 is
the step-boundary audit (`plan-authoring` §6b/§6c). §7 is the mutation battery each gate
owes. §8 is the refuted / unworkable list. §9 the residual gaps and the open questions for
the user.

---

## 0. Headline findings (read these before designing anything else)

| # | finding | consequence |
|---|---|---|
| **F1** | ⭐⭐ **A byte-identity wall for CS3 already exists and is CHEAP.** `tests/sn/solve/test_affine_carve_bit_identity.py` under `-W error::tests.sn.regression._regression_assert.DriftWarning` is a 1-ULP gate on three drivers (2-D windowed SI, 2-D Krylov, 1-D slab SI). `[M]` **3 passed in 1.60 s**; positive control (1-ULP perturbation of the loaded baseline) **3 failed in 1.65 s**. | Do NOT mint a new bit-identity instrument for those three. Escalate this one. §1.2 |
| **F2** | ⭐⭐ **"SI computes the diagnostic trajectory from flat norms" is AMBIGUOUS, and one reading silently moves ρ by 0.47 %.** Today ρ is a ratio of `Displacement.l2` = `space.norm(leaf.values)` on the **INTERIOR LEAF ONLY**. `[M]` interior-metric vs interior-flat: max rel **2.29e-16** (≤1 ULP). Interior-metric vs **whole-composite** flat (`_l2_norm(displacement)`, the spelling the SI loop already has in hand): max rel **4.71e-3**. | The relocation must state WHICH norm and WHICH block. §2 pins it; the gate carries the 4.71e-3 control. |
| **F3** | ⛔ **The DSA input-type refusal has NO negative test.** `[M]` `grep -rn "DSACorrection" tests/` → 11 hits, **zero** `pytest.raises` against `DSACorrection.apply`. Step 2 rewrites exactly that `isinstance(interior, (AngularDisplacement, AngularFlux))` arm. | The replacement's teeth are NET-NEW, not migrated. §4.5 writes the test the guard never had. |
| **F4** | ⛔ **DSA-accelerated SI has no old-vs-new VALUE gate.** Its only invariance gate is `test_dsa_acceleration.py`'s `rtol=1e-7` FP-invariance against plain SI — a *proxy* in exactly the `vv` #12 sense (it compares two runs of the *same* carved code). | §1.3 specs a 4th case for the `test_affine_carve_bit_identity` harness. |
| **F5** | ⛔ **The adjoint solve routes typed carriers through `KEigenvalue.measure_stopping_criteria`, which evaluates `flux_distribution - flux_old` — i.e. it MINTS A DISPLACEMENT today** (`iteration.py:1526`, carrier = `FullField`/`CoupledField` per `solver.py:2930-2935`). It has no stored-value gate either; its gates are THEOREM rows (`k_adj == k_fwd`, `atol=1e-8`). | §1.4. The SN *forward* eigen path is NOT exposed the same way — `SNSolver.measure_stopping_criteria` takes bare `np.ndarray` (`solver.py:1769`). |
| **F6** | ⭐ **The fiber refusal really does survive the flip — measured, not argued.** Simulating the post-carve algebra by calling `Field.__add__` directly on two cross-mesh same-shape fluxes: `[M]` REFUSED with `ValueError: <Leaf> arithmetic across distinct SNMesh instances is forbidden — the field is mesh-bound` on `AngularFlux`, `ScalarFlux`, `AngularBoundaryFlux`, `HarmonicMomentFlux`; same-mesh add returns the leaf type. ⚠ The **space** gate does NOT catch it: `[M]` `a.space == b.space` is `True` (`FunctionSpace.__eq__` is `(name, shape)`) while `a.space is b.space` is `False`. Only the mesh arm refuses. | §4.3 pins it with the message fragment `across distinct SNMesh instances`. |
| **F7** | ⭐⭐ **A production-tier cone violation exists and is cheap** — step 4 does NOT need the cell-level witness the brief names. `[M]` `solve_sn_fixed_source` on a 4-cell / 40 cm / Σ_t = 10 slab (Δx·Σ_t = 100, S2, unit source in the FIRST cell only) converges to `min ψ = −6.399383e-01`, **2 of 8 entries negative**, `min φ = −8.438399e-01`. A benign sibling (nx = 2, same materials) gives `min ψ = +2.181405e-01`. | §5 uses the public entry for both legs. The brief's `TestPositivityFailure` is a CELL-level witness (`strat.update(...)` → a bare ndarray), not a field — §8.2. |
| **F8** | ⛔ **CS2 will legitimately RED the ρ pin.** `[M]` `angular_flux`'s space has `inner_product_weights is None` today (Euclidean), which is *why* F2's metric-vs-flat gap is ≤1 ULP. Recomputing the same trajectory under a physical `V_cell × w_n` metric moves ρ by **1.12e-3** relative — 9 orders above the pin's `rtol=1e-12`. | The pin is a CS3-scoped instrument and its docstring says so. CS2's F3 metric relocation owns re-deriving it (or the diagnostic must be *defined* on the Euclidean norm). §2.4 |

---

## 1. Value-neutrality harness (`vv` #12: no proxies — a direct old-vs-new VALUE comparison per consumer)

### 1.0 Why the carve is expected byte-identical, and why that is not evidence

The two arithmetic expressions are *character-identical* before and after:

| op | today (`FluxRole`) | after (base `Field`) |
|---|---|---|
| `ψ − ψ'` | `_check_partner(other)` then `_mint_displacement(self.values - other.values)` | `_check_partner(other)` then `replace(self, values=self.values - other.values)` |
| `ψ + Δ` | `_check_torsor_partner(other)` then `replace(self, values=self.values + other.values)` | `_check_partner(other)` then `replace(self, values=self.values + other.values)` |

So the SI update `psi + corrector.apply(psi - psi_prev)` (`iteration.py:780`) evaluates the
identical numpy ops on the identical buffers; only the RESULT TYPE moves. Composite algebra
is fully member-delegating (`full_field.py:384-390` — `__add__`/`__sub__` are
`_map_binary(other, lambda a, b: a + b)`), so leaf-level neutrality propagates.

⛔ **That argument is a hypothesis, not a gate.** `vv` #12 exists because a neutrality claim
holds only for the ONE contract it was proven against. The table below is the proof
obligation.

### 1.1 Consumer × gate matrix

| # | production consumer | flux-algebra exposure | existing VALUE gate | verdict |
|---|---|---|---|---|
| C1 | **SI, 1-D slab** (`AngularFlux` bulk ⟹ `AngularDisplacement`) | `iteration.py:787` mint; `:780` corrector torsor add | `test_affine_carve_bit_identity.py::[si_slab_2g_het]` — stored `.npy`, `SAFETY × conv_tol = 1e-11`, `DriftWarning` on ANY movement | ✅ **escalate** (F1) |
| C2 | **SI, 2-D windowed** (`HarmonicMomentFlux` bulk ⟹ `MomentDisplacement`) | same, moment leaf | `…::[si_2d_p1_aniso_het]` | ✅ **escalate** |
| C3 | **Krylov** | none on the flat path (roles erased at `to_flat`/`from_flat`); exposure is only the template rebuild | `…::[krylov_2d_p1_aniso_het]` | ✅ **escalate** — and it is the *control*: a red here means the carve leaked into a path that has no flux dunder in it |
| C4 | **DSA-accelerated SI** | `dsa.py:635` isinstance dispatch; `:680-696` mints `AngularDisplacement` + `AngularBoundaryDisplacement`; consumed by `iteration.py:780` | ⛔ **NONE.** `test_dsa_acceleration.py` compares DSA against plain SI at `rtol=1e-7` — both sides carved | ⛔ **new case, §1.3** |
| C5 | **`power_iteration`, SN forward keff** | outer is bare-ndarray (`solver.py:1769-1793`, `np.linalg.norm(flux_distribution - flux_old)`); exposure is the INNER SI only | `tests/sn/regression/test_dd_regression.py` — 13 cases. `[M]` **11 of 13 are bit-exact at HEAD**; the 2 pre-drifting rows are exactly `cyl_1g_homogeneous_folded_2x4_dd_n20` and `…_4x8_…` (`scalar_flux` 40 679 ULP / 5.32e-12 and 549 721 ULP / 7.19e-11) | ✅ **escalate with those 2 deselected**, §1.5 |
| C6 | **Adjoint solve** (`solve_sn_adjoint`) | `KEigenvalue.measure_stopping_criteria` `flux_distribution - flux_old` on a TYPED carrier ⟹ mints a displacement composite | ⛔ **no stored-value gate.** `test_sn_adjoint_certification.py` is THEOREM-class (`k_adj == k_fwd`, `atol=1e-8`) | ⛔ **new case, §1.4** |
| C7 | **Coupled / curvilinear SI** (`CoupledField`, the 2 RC leaves) | `_flux_displacement_leaf` recursion into `systems[0]`; RC displacement leaves minted per block | wiring: `test_psi_half_coupling.py::test_g_d1_7_…` (diagnostic non-empty). values: the `sphere_*` / `cyl_*` DD snapshots | ⚠ wiring gate **re-points** (§3); values covered by C5's suite |

### 1.2 The bit-identity wall — the command, and its measured positive control

```
.venv/bin/python -O -m pytest tests/sn/solve/test_affine_carve_bit_identity.py -q \
  -W error::tests.sn.regression._regression_assert.DriftWarning
```

`[M]` at `000cf144`: **3 passed, 1.60 s**.

Positive control (`vv` #17 — a clean reading carries no information without one). A pytest
plugin that monkeypatches `np.load` to advance the FIRST element of every loaded baseline by
one ULP, and `raise`s at `sessionfinish` if it perturbed nothing (the harness asserts its own
installation):

- escalated: `[M]` **3 failed, 1.65 s** — all three cases.
- unescalated: `[M]` **3 passed, 7 warnings**, e.g.
  `si_slab_2g_het:scalar_flux: scalar_flux drifted 1 ULP / 1.90e-16 rel (within tol 1.0e-11)`.

⟹ the escalation is a genuine 1-ULP instrument, the `-W` string parses (`vv` Mode-8 EIGHTH
class: gate the STRING, not the API), and the unescalated run is a *magnitude reporter* rather
than a pass/fail — both readings are useful and they are different claims.

### 1.3 The missing DSA case (C4)

**Add a fourth case to the EXISTING harness** — never a second capture mechanism (that file
already owns `--capture-baseline`, `_baseline_path`, `_capture_or_assert`, and the
regeneration-history discipline; a parallel copy is a Pattern-2 twin).

- `case = "dsa_slab_2g_het"`: `_build_slab()` unchanged (2 G, fuel|moderator, P1-aniso,
  vacuum), `inner_solver="source_iteration"`, `acceleration="dsa"`, same `_INNER_TOL` /
  `_MAX_INNER`.
- ⚠ **Capture BEFORE step 1 lands**, at `000cf144`, with `--capture-baseline`; commit the
  two `.npy` files in the same commit as the case.
- ⚠ Prefer a **reflective** wall on one side (`bc_left=BC("reflective")`) or add a second
  case with one — `dsa.py`'s trace arm is documented "load-bearing under the lagged
  reflective gain; inert on vacuum faces", and a vacuum-only case makes the whole
  `AngularBoundaryDisplacement` mint (the block whose TYPE step 2 changes) a
  provable non-catcher (`vv` #20 row-inflation: a row that cannot see the varied thing).
  **This is the activation requirement for C4** — without it the DSA case pins the bulk arm
  only and the boundary re-type ships ungated.

### 1.4 The missing adjoint case (C6)

`solve_sn_adjoint` on a small heterogeneous ≥2 G slab, converged `k_adj` + `psi_star`
System-A bulk values stored the same way. Two reasons it cannot be skipped:

1. It is the ONLY consumer where a typed-carrier difference is taken **outside**
   `SourceIteration` (`iteration.py:1526`), so a relocation that fixes up the SI loop and
   forgets `KEigenvalue` shows here and nowhere else.
2. `solver.py:2960` carries a comment justifying an OMITTED balance-defect clause on the
   grounds that "the `1/k` scaling of a field crosses the affine-torsor arithmetic rules
   (a flux state is not a vector)". Memo D §2 already records this as false *today* (scalar
   scaling is legal); after CS3 it is doubly void. **That comment and issue #353 are in the
   carve's doc blast radius** (`coding-standards`: a present-tense-false justification is a
   MUST-FIX; `vv` anti-#21: grep a WINDOW, the subject and its negation wrap).

### 1.5 The broad value wall (C5, C7)

```
.venv/bin/python -O -m pytest tests/sn/regression/test_dd_regression.py -q \
  -W error::tests.sn.regression._regression_assert.DriftWarning \
  --deselect "tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_1g_homogeneous_folded_2x4_dd_n20]" \
  --deselect "tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_1g_homogeneous_folded_4x8_dd_n20]"
```

`[M]` at `000cf144`: **11 passed, 52.7 s** (the un-deselected run is `2 failed, 11 passed`,
and those 2 are the characterised pre-drift, NOT a CS3 signal). Characterising rather than
counting is what makes every post-carve red attributable with no triage.

⚠ Do NOT re-baseline anything in this suite as part of CS3. If a row moves, the carve is not
neutral and the finding is the point.

---

## 2. The ρ-trajectory equivalence gate — ⏹ WRITTEN, GREEN, MUTATION-VERIFIED

**File**: `tests/numerics/test_si_diagnostic_trajectory.py` (new, 5 tests).
**Run**: `.venv/bin/python -O -m pytest tests/numerics/test_si_diagnostic_trajectory.py -q`
→ `[M]` **5 passed in 0.46 s**. `pyright` → **0 errors, 0 warnings**.

### 2.1 The fixture and its configuration (`plan-authoring` §4)

2 groups · **heterogeneous** (2 zones of 20 cells) · asymmetric downscatter ·
`c = Σ_s/Σ_t = 0.99` in **both groups of both zones** · `Σ_t = (1.0, 2.0)` fuel /
`(1.5, 1.0)` moderator · slab `[0, 10] cm` (10 mfp) · `nx = 40` · Gauss-Legendre
`n_ord = 8` · vacuum|vacuum · uniform isotropic unit source · `max_iter = 12`,
`tol = 1e-14` (so the stop never fires ⟹ 11 ratios, deterministically).

⚠ **2G heterogeneous is not decoration.** A 1-group homogeneous slab makes a
per-group norm and a total norm the same number, so a relocation reducing over the
wrong axis would be Mode-12 invisible. The sibling ρ≈c reference gate is 1G, and
legitimately so (a rate claim is flux-shape-independent) — this pin is not.

### 2.2 What is frozen

`contraction_ratios` (11 values, `[M]` the trajectory CROSSES 1.0 in the transient —
`0.97593, 1.03588, 1.02699, …, 0.97127` — so it is a fingerprint, not a constant),
`‖Δψ‖ = 8.848958875813311`, `‖Δψ‖/(1−ρ) = 308.0488097665712`,
`where_largest(3) = [(3,1,24), (4,1,26), (4,1,25)]`.

### 2.3 The tolerance, derived rather than chosen

`rtol = 1e-12`, `atol = 0`. Between two measured floors:

| what | `[M]` max relative | verdict |
|---|---|---|
| swap `space.norm` → `np.linalg.norm` on the interior leaf | **2.29e-16** (0–1 ULP/step) | below the pin — a legitimate, DECLARED blindness |
| use the WHOLE composite's flat norm (interior ⊕ boundary trace) | **4.71e-3** | 9 orders above the pin — the error it exists to catch |

### 2.4 Battery (in-process plugins; the harness `raise`s if it perturbs nothing)

| # | mutation | `[M]` result |
|---|---|---|
| M1 | ρ from the whole composite's flat norm | **2 failed** — ρ pin AND ‖Δψ‖ pin |
| M2 | `Field.l2` → `np.linalg.norm(values)` | **5 passed** — the declared ≤1 ULP blindness, measured |
| M3 | one spurious leading ratio (cadence drift) | **1 failed** — the length leg |
| M4 | `where_largest` ravels first (layout loss) | **1 failed** — the map leg only |
| M5 | every norm × (1+1e-9) — positive control | **1 failed** — the ‖Δψ‖ leg **ONLY** |

⛔ **M5 is a stabiliser finding about the gate itself, and it is recorded in the
module docstring**: the ρ trajectory is INVARIANT under any factor applied
uniformly to the norm (ρ is a ratio; a common scale cancels exactly — `vv` Mode
12). That whole class is caught only by the separate ‖Δψ‖ pin, which is why the
two must never be folded together.

### 2.5 The measured CS2 hazard, and the decision it forces

`[M]` recomputing the same trajectory under a physical `V_cell × w_n` metric moves
ρ by up to **1.12e-3** relative. So **CS2 will legitimately RED this gate**. CS3
must therefore RULE, in the relocation itself, which norm the diagnostic is defined
on:

- **(a) space norm** — the diagnostic follows the physics metric; CS2 owns
  re-deriving these numbers, with a regeneration note in the style of
  `test_affine_carve_bit_identity.py`'s.
- **(b) Euclidean/flat norm** — metric-independent by construction, these numbers
  are permanent, and the diagnostic stops being a Hilbert-space quantity.

Either is defensible; leaving it undecided is not, because the gate reads as a
regression under (a) and as a contract under (b). Recorded as **Q1** in §9.

---

## 3. Test-migration map

Verdict vocabulary: **RE-SPELL** (same claim, new syntax) · **RE-POINT** (same claim,
new surface) · **RETIRE** (the claim ceases to exist) · **STRENGTHEN** (the claim gets
better because the carve removes a workaround).

### 3.1 The foundation battery — `tests/numerics/test_affine_flux_algebra.py` (10 tests, `foundation`)

⚠ **This one file is cut by the step-1/step-2 boundary**: 3 of its 10 tests are
step-1 material (the diagnostics) and 7 are step-2 (the algebra). Landing step 1
alone does not break it — but see §6.1, because leaving `Displacement`'s three
methods alive while the iteration layer grows its own copy is a transient
Pattern-2 twin, and these 3 tests would then pin the DEAD copy.

| test | step | verdict |
|---|---|---|
| `test_mint_stores_raw_difference_and_torsor_round_trips` | 2 | **RE-SPELL** → `ψ₂ − ψ₁` returns the FLUX type; the round-trip `ψ₁ + (ψ₂ − ψ₁) ≈ ψ₂` survives verbatim (still `nulp=8` — `a + (b − a) ≠ b` bit-for-bit). Keep the bit-exact "stores the raw difference" leg; INVERT the type assertion (was: MUST NOT be the flux type). |
| `test_displacement_telescoping` | 2 | **RE-SPELL** — path-independence of V is now a statement about the SAME type. Unchanged arithmetic. |
| `test_displacement_scalar_action_and_zero` | 2 | **RE-SPELL**; add the NEW law `ψ + 0 == ψ` where `0` is a *flux*, which was unspellable. |
| `test_flux_add_flux_forbidden_but_torsor_allowed` | 2 | **RETIRE the negative leg, INVERT into §4.1** — `flux + flux` becomes the headline legal op. |
| `test_displacement_add_displacement_returns_displacement` | 2 | **RETIRE** (collapses into the re-spelled telescoping row). |
| `test_affine_combination_partition_of_unity` | 2 | **RETIRE** — `affine_combination` retires with 0 production callers. ⚠ Its *content* (`ω·ψ_new + (1−ω)·ψ_old` equals the relaxation form) is a real claim about relaxation; keep ONE row of it re-spelled as plain algebra in §4.1, or the carve silently drops a blend law. |
| `test_subtraction_typed_and_guards_intact` | 2 | **RE-POINT** — the cross-class and cross-mesh legs are the fiber discipline and MUST survive; the "typed displacement" leg inverts. This is §4.3's ancestor. |
| `test_contraction_ratio_of_scaled_displacement` | 1 | **RE-POINT** to the iteration layer. Keep: a geometric step `Δ⁽ⁱ⁺¹⁾ = ρ·Δ⁽ⁱ⁾` must recover ρ exactly. |
| `test_true_error_estimate_amplifies` | 1 | **RE-POINT**. Keep BOTH legs (`ρ=0.9 ⟹ 10×`, and `ρ≥1 ⟹ ValueError`) — the refusal is the `vv` #11 negative. |
| `test_where_largest_locates_the_peak` | 1 | **RE-POINT**. |

The module's own docstring, header maths, and its reference to
`issue_208_affine_algebra_verification_spec.md` all become historical — keep them
past-tense (`coding-standards`: past-tense history stays, present-tense falsehood is
the MUST-FIX). Rename the file to something the cone ontology owns.

### 3.2 The 7 files pinning `flux + flux` (all step 2)

| file | site(s) | verdict |
|---|---|---|
| `tests/numerics/test_affine_flux_algebra.py` | §3.1 | as above |
| `tests/transport/fields/test_angular_flux.py:106-114` | `raises(TypeError, match="affine_combination")` + torsor | **RE-SPELL → §4.1**; its `:155` comment about where the cross-mesh guard lives is the §4.3 claim — keep and re-point |
| `tests/transport/fields/test_angular_boundary_flux.py:98-113, 146-160, 287, 330` | mint + distributive + per-face difference | **RE-SPELL** (the per-face-difference and distributive rows survive unchanged in content) |
| `tests/sn/primitives/test_harmonic_moment_flux.py:226-245, 350-380` | the `iso + aniso` recombination refusal + the moment torsor | **RE-SPELL — and re-examine the CLAIM.** `:226` says the ℓ-disjoint recombination `iso + aniso` is "FORBIDDEN by the #208 affine gate". After the carve it is LEGAL and it is exactly the right spelling. That is a capability the carve unlocks; the test should assert it, not its refusal. |
| `tests/transport/test_timed_full_field.py:176-186, 230-240, 258` | composite propagation | **RE-SPELL**; keep the composite-level distributive row |
| `tests/sn/primitives/test_typed_source_sinks.py:432-461` | the CONTRAST row (source-sink `+` is closed; flux `+` is not) | ⚠ **RE-POSE, do not delete.** Its content is "the rate-density family and the state family have DIFFERENT algebra". After CS3 they have the same additive algebra, so the contrast collapses — the surviving distinction is the cross-class gate (`AngularFlux + AngularSourceSink` still refuses). Re-point at that. |
| `tests/sn/sweep/cartesian_2d/test_2d_l2_matvec_correctness.py:283-340` | a hand-built `.values`-level workaround **because** `α·u + β·v` was illegal | ⭐ **STRENGTHEN.** Delete the workaround and the `raises(TypeError)`; the test becomes the direct typed linearity statement `A.apply(α·u + β·v) == α·A.apply(u) + β·A.apply(v)`. This is the carve's best advertisement, and it is a gate that gets STRONGER (it moves from raw arrays to the typed algebra). |

`tests/sn/sweep/core/test_ordinate_scan.py` is a **FALSE POSITIVE** of the concept
grep: its "affine" is the DD scan recurrence being affine in `(b, ψ₀)`, nothing to
do with the field algebra. Do not touch it. (`coding-standards`: a concept grep's
hits must be triaged by MEANING before any is called a site.)

### 3.3 The ~16 affine raise-sites — the exact denominator

`[M]` `grep -rn "raises(TypeError" tests/ -A2 | grep -c "affine\|origin"` → **16**, in
10 files: `test_affine_flux_algebra` 6 · `test_harmonic_moment_flux` 2 · one each in
`test_timed_full_field`, `test_radial_characteristic_field`, `test_composite`,
`test_scalar_boundary_flux`, `test_angular_flux`, `test_angular_boundary_flux`,
`test_typed_source_sinks`, `test_radial_characteristic_split_leaves`. All step 2.
The three not already covered above:

- `tests/transport/test_composite.py` — container-level delegation of the gate.
  **RE-SPELL** (the delegation claim survives; only what is delegated changes).
- `tests/transport/test_radial_characteristic_field.py:172-174` and
  `tests/sn/mesh/test_radial_characteristic_split_leaves.py:150-166` — assert the
  composite mints RC interior/boundary displacements and that the two leaves are
  DISTINCT classes. **RETIRE the type rows** (their subject is the leaf pair) and
  keep the "the composite propagates blockwise" claim re-spelled.
- `tests/transport/fields/test_scalar_boundary_flux.py:199-207` — **RE-SPELL**.
  ⚠ Build note: `ScalarBoundaryFlux.zeros_on` needs a `DiffusionMesh`, not a bare
  `SNMesh` (`[M]` `_bases.py:1117` refuses: "mesh carries no scalar trace").

### 3.4 The 9 Displacement-referencing files, and two special cases

`[M]` `grep -rln "Displacement" tests/` → 9 files. Six are covered above. The rest:

| file | verdict |
|---|---|
| `tests/sn/solve/test_flux_displacement_diagnostics.py` (4 tests, `l1`) | ⭐ **RE-POINT, never RETIRE — this is the ρ REFERENCE anchor** (ρ ≈ c, Adams & Larsen 2002) that keeps §2's RECORD pin honest. Its `_asymptotic_rho` reads `si.contraction_ratios` and `si.last_displacement`; both re-point with step 1. Its 1-group fixture stays legitimate (a rate claim is flux-shape-independent). |
| `tests/sn/operators/test_psi_half_coupling.py:3613-3636` (`test_g_d1_7_…`) | ⭐ **RE-POINT — and it is a load-bearing catcher.** It pins that the diagnostics record on a COUPLED (sphere, `CoupledField`) iterate. Today `_flux_displacement_leaf` duck-types on `contraction_ratio`; after step 1 no leaf has that method, so a naive relocation returns `None` and the diagnostics go SILENT on the coupled path — exactly this test's stated failure mode. Its replacement must be a STRUCTURAL walk (`.interior`, else `.systems[0].interior`). |
| `tests/sn/acceleration/test_dsa_low_order.py:259-277` (`test_d8_displacement_reduction_is_the_tangent_map`) | **RETIRE.** After the carve it is a character-for-character duplicate of `test_d8_restriction_is_the_frame_moment_row` two functions above (flux − flux → flux ⟹ `integrate_angular` → `ScalarFlux` ⟹ the same body on the same values). `[M]` it carries NO `verifies` marker of its own (the two `sn-dsa-restriction` markers sit on lines 227 and 245), so no marker migration is owed here. |

### 3.5 Marker migration (`coding-standards`: retirement = marker migration)

`[M]` grep of the 12 affected files: **zero `catches("ERR-NNN")` markers** anywhere in
the CS3 blast radius, and exactly **two `verifies(...)`** — both
`sn-dsa-restriction` / `sn-dsa-consistent-low-order` in `test_dsa_low_order.py`, and
both on tests that SURVIVE. ⟹ **marker migration for CS3 is empty, and that is a
measured fact rather than an omission.** Module-level marks that must be carried onto
any renamed/replacement module: `foundation` (7 files), `l1` (the diagnostics
anchor), `l2` (the DSA acceleration gates).

⚠ Separately noted, not CS3's job: `tests/transport/fields/test_angular_flux.py`,
`…/test_angular_boundary_flux.py`, `…/test_scalar_boundary_flux.py` and
`tests/transport/test_timed_full_field.py` carry **no V&V marker at all** (`[M]` no
`pytestmark`, no `@pytest.mark.*`). A file being edited anyway is the cheap moment
to tag them `foundation`.

### 3.6 Doc surface — and the grep trap that will bite

5 source pages spell `torsor`: `field_algebra.rst` (the 602-line page whose
§"Why affine" at `:242` is the OVERTURNED argument), `operator_algebra.rst` (11 ref
sites), `methods/sn/acceleration.rst:887`, `methods/sn/curvilinear_one_group.rst:4170`,
`verification/matrix.rst:875`.

Four labelled equations are APIs (`coding-standards`: a labelled equation is an API
the moment anything writes `:eq:`): `affine-torsor-algebra` (`:176`),
`affine-contraction-ratio` (`:399`), `affine-true-error` (`:418`),
`affine-typed-residual-eq` (`:476`). `[M]` external `:eq:` citers of the four:
`operator_algebra.rst:3402` and `:3955` only.

⛔⛔ **THE TRAP: `affine-bc-form` is a DIFFERENT label and it must not be touched.**
It is the affine BOUNDARY law `γ₋ψ = R G γ₊ψ + q`, unrelated to the flux torsor, and
`[M]` it has ~**18** `:eq:` citers across `foundations/boundary_conditions.rst`,
`methods/sn/boundary_conditions.rst` and `verification/sn.rst`. A sweep grepping
`affine` will return them all. **Grep `torsor` and the four exact label names; never
the bare word `affine`.**

Also in the doc/comment blast radius, all present-tense-false after the carve:
`orpheus/numerics/coupled_system.py` (`:108`, `:239`, `:292`), `numerics/vector.py:55-58`,
`transport/__init__.py:31-32`, `transport/source_sinks/harmonic_moment_source_sink.py:27-42`
(its whole rationale is "the role split IS the design"), `sn/solver.py:2960` + issue
**#353**, and — outside the tree — the `coding-elegance` skill's anti-pattern **#18**,
which cites `FluxDisplacement` as a POSITIVE precedent and would otherwise teach the
retired ontology.

---

## 4. New-algebra gates — what pins the laws the carve CREATES

Every row declares its claim layer, its pillar, and the mutation that reddens it
(`AGENT.md` §0.5 / §1.5). All are `foundation` (software invariants, no `:label:`),
all use `np.testing.assert_*` or explicit `raise` (`-O`).

### 4.1 `flux + flux` is LEGAL and arithmetically correct — the headline

**Claim layer**: none (a software algebra law). **Pillar**: closed-form — numpy `+`
on the same buffers, structurally independent of the field layer.
**Home**: the re-pointed foundation battery, parameterized over ALL 7 flux leaves
(`AngularFlux`, `ScalarFlux`, `HarmonicMomentFlux`, `AngularBoundaryFlux`,
`ScalarBoundaryFlux`, `RadialCharacteristicInteriorFlux`,
`RadialCharacteristicBoundaryFlux`).

| leg | assertion |
|---|---|
| a | `(ψ₁ + ψ₂).values` is `array_equal` to `ψ₁.values + ψ₂.values` — bit-exact, single numpy add |
| b | `type(ψ₁ + ψ₂) is type(ψ₁)` — the sum stays in the SAME leaf type (this is the whole ontology change) |
| c | commutativity `ψ₁ + ψ₂ == ψ₂ + ψ₁` and associativity to a few ULP |
| d | `ψ + 0` is `array_equal` to `ψ` where `0` is a *flux* built by `zeros_on` — **the origin exists**, the exact statement the "no natural zero" doctrine denied |
| e | the retired blend: `0.7·ψ₂ + 0.3·ψ₁` equals the old `affine_combination` value — carries `affine_combination`'s only real content forward without its Σλ=1 ceremony |

⚠ **Leg (d) is the one to write carefully.** "A zero flux constructs" is already true
today (`zeros_on` is used everywhere); what is NEW is that it is an ADDITIVE
IDENTITY. Assert `ψ + zero is not ψ` **and** `array_equal(values)` — otherwise a
`return self` short-circuit would pass while breaking the copy contract.

**Mutation**: replace `Field.__add__`'s body with `self.values - other.values`.
Expect: legs (a), (c), (d), (e) red; leg (b) stays green (the TYPE is unchanged) —
which is exactly why (b) alone is not the gate.

### 4.2 `flux − flux` returns the SAME type and is SIGNED

| leg | assertion |
|---|---|
| a | `type(ψ₁ − ψ₂) is type(ψ₁)` — **the inverted assertion**; today's battery raises if this holds |
| b | values `array_equal` to `ψ₁.values − ψ₂.values`, and the result has **negative entries** (build ψ₁ < ψ₂ somewhere) — the SIGNED claim, and it is the cone's own boundary: a difference of cone members leaves the cone |
| c | telescoping `(ψ₂−ψ₁) + (ψ₃−ψ₂) ≈ ψ₃−ψ₁` to `nulp=8` |
| d | `Displacement` is **not importable** (post-step-3): a `pytest.raises(ImportError)` / `importlib.util.find_spec("orpheus.transport.displacements") is None` row is the retirement's own gate |

⭐ Leg (b) is the seam between §4 and §5: it states in one assertion that V is signed
and that K is a strict subset — the exact reason cone membership must be a predicate
and never a constructor invariant.

### 4.3 The fiber refusal is RETAINED — the charter's owed negative control

`[M]` **already proven runnable** (F6): post-carve `Field.__add__` on two cross-mesh
same-shape fluxes REFUSES.

| leg | assertion |
|---|---|
| a | **NEGATIVE** — `ψ_meshA + ψ_meshB` raises `ValueError` with `match="across distinct SNMesh instances"`. Parameterize over all 7 leaves. |
| b | **POSITIVE** — `ψ_meshA + ψ'_meshA` returns the leaf type (`vv` #11: a guard that refuses everything also passes a negative-only test) |
| c | **THE DISCRIMINATOR** — assert in-test that `meshA is not meshB` **and** `ψ_meshA.space == ψ_meshB.space` is `True`. This pins that the refusal comes from the MESH arm and not from the space gate; without it the row silently degrades the day `FunctionSpace.__eq__` starts distinguishing instances. |
| d | cross-class still refuses: `AngularFlux + AngularSourceSink` → `TypeError` (Layer 1 survives; this is what `test_typed_source_sinks.py`'s contrast row re-poses onto, §3.2) |

**Mutation**: delete the mesh arm in `_bases.py` (`self.mesh is not other.mesh`).
Expect leg (a) red on every leaf, legs (b)–(d) green. ⚠ Run this mutation — the
whole "the fiber discipline survives" claim rests on ONE `is` comparison in a base
class, and `_check_torsor_partner`'s own mesh arm (`_flux_role.py:196`) retires with
the mixin, so the surviving copy is not the one CS3 readers will be looking at.

### 4.4 Scalar algebra unchanged (the regression floor)

`c·ψ`, `ψ/k`, `−ψ` return the leaf type with `array_equal` values. Cheap, and it is
the row that catches a carve that "simplified" `FluxRole` away by dropping the whole
mixin file including something else. Also assert `ψ/k` works for the eigenvalue
normalisation — the operation `sn/solver.py:2960` wrongly claims is forbidden.

### 4.5 The DSA input-type refusal — teeth that are NET-NEW (F3)

`[M]` no `pytest.raises` anywhere targets `DSACorrection.apply`. So its guard has
never had a negative test, and step 2 rewrites it.

| leg | assertion |
|---|---|
| a | **NEGATIVE** — `DSACorrection.apply(<composite with a HarmonicMomentFlux interior>)` raises `TypeError`, `match=` a SHORT distinctive fragment of the shipped message (`"moment-windowed"` — never the full sentence; `coding-standards`: consumers pin substrings) |
| b | **POSITIVE** — the same corrector applied to a full-angular composite returns a composite whose interior is the flux type and whose boundary block is non-zero on a REFLECTIVE fixture |
| c | leg (b) on a VACUUM fixture as the declared control: the trace arm is documented inert there, so assert the boundary block is (near) zero and say in the docstring that this row cannot see the trace arm |

⚠ This must be written BEFORE step 2 touches the isinstance arm, and it must be RED
against a mutation that simply deletes the guard.

### 4.6 The `#331` closure — operator linearity through the typed algebra

`[M]` at HEAD the two gain leaves refuse a displacement, verbatim:

- `ScatteringOperator.apply: unsupported input type AngularDisplacement; expected TimedFullField, ScalarFlux, Ang…`
- `SNBoundaryOperator: the input composite's boundary must be an AngularBoundaryFlux trace; got AngularBoundaryDi…`

After CS3 both are unspellable. The done-when is a POSITIVE gate, and it is a claim
that could not previously be written at all:

```
S.apply(ψ₁ + ψ₂)  ==  S.apply(ψ₁) + S.apply(ψ₂)      # additivity, typed
S.apply(c · ψ)    ==  c · S.apply(ψ)                  # homogeneity, typed
S.apply(ψ₁ − ψ₂)  ==  S.apply(ψ₁) − S.apply(ψ₂)      # the #331 case, typed
```

⚠ **Activation leg, mandatory** (lessons L40c): measure `‖S.apply(ψ)‖ > 0` for the
random ψ and assert it. A linearity row on an operator that is the zero morphism on
the chosen fixture holds with both sides structurally zero and no input can red it.
Apply the same to `B` — on a fixture where BOTH faces are `PrescribedInflow`, `B` IS
the zero morphism (P3 collapsed it), so the fixture must be prescribed-on-one-face /
REFLECTIVE-on-the-other.

`tests/sn/operators/test_declared_law_is_linear.py` then **RE-SPELLS**: its
base-point trick (`B(ψ₁+σ) − B(ψ₂+σ) == B(ψ₁) − B(ψ₂)` with σ a displacement) exists
only because plain additivity was unspellable. Its assertions are value-level, so the
re-spell is mechanical — but keep the module header's explanation as past-tense
history, because it records WHY the trick existed.

---

## 5. Cone-predicate gate spec (step 4)

### 5.1 What the predicate is, and what the test does NOT claim

The predicate observes; **nothing enforces**. There is no fixup, no clipping, no
warning anywhere in the SN spatial path (`[M]` memo D §3.2, re-checked: `grep -rn
"fixup\|clip" orpheus/sn/ orpheus/transport/spatial/` returns prose only). So the
gate's claim is exactly:

> *the predicate correctly reports whether a given field lies in K*

and explicitly **NOT**: that production keeps fields in K, that a violation is
handled, that DD is repaired, or that a violating solve is refused. The docstring
must say all four negatives, or the next audit will read the green as "cone
preserved".

### 5.2 The precedent, read rather than cited (`plan-authoring` §1)

`CrossSectionField`'s cone-as-predicate doctrine (`cross_section_field.py:30-43`) is
the right model and its battery `TestCrossSectionConeAlgebra`
(`tests/transport/fields/test_coefficient_fields.py:75-120`, 5 tests) is the right
SHAPE: closed under `+`, closed under `λ≥0·`, has an origin, is a vector space not a
torsor, units.

⛔ **Correction to "model it on X": there is NO predicate to copy.** `[M]` `grep -rn
"def is_nonnegative\|def in_cone\|def is_in_cone" orpheus/transport/fields/
orpheus/numerics/field.py` → **0 hits**. The precedent supplies the DOCTRINE and the
test shape; the method itself is net-new. (Also note the precedent battery uses bare
`assert` — legal there because pytest rewrites collected modules, but the new battery
should use `np.testing.assert_*` / `raise` to match house style.)

### 5.3 The two legs, both through the PUBLIC entry (F7)

| leg | fixture | `[M]` at HEAD | assertion |
|---|---|---|---|
| **positive** | `solve_sn_fixed_source`, slab `nx=2`, `width=20`, `Σ_t=10`, `c=0.5`, S2, source in cell 0 | `min ψ = +2.181405e-01` | the predicate accepts the converged `angular_flux` AND the `scalar_flux` |
| **negative** | identical materials/quadrature/BC, `nx=4`, `width=40` (Δx·Σ_t = 100) | `min ψ = −6.399383e-01`, **2 of 8 entries negative**; `min φ = −8.438399e-01` | the predicate REJECTS, and reports WHERE (`vv` anti-#14: return the structure, not a `bool` — a predicate that returns the offending indices makes its own correctness assertable) |

⭐ Why this pair is strong: the two rows differ in **one** parameter (`nx`), so the
negative leg cannot be explained by a materials or quadrature difference, and both
run the same production entry. `[M]` both solves converge.

### 5.4 The unit-level legs (no solver)

- `zeros_on(mesh)` IS in K (the cone contains its apex) — and after §4.1(d) the same
  object is the additive identity, so K's apex and V's origin coincide. State it.
- K is closed under `+` and under `λ ≥ 0 ·`; **not** closed under `−` (§4.2 leg (b)'s
  signed difference is rejected by the predicate) and **not** under `λ < 0 ·`.
- a hand-built field with exactly ONE negative entry at a known index is rejected AND
  the reported index is that one — the sharpest cheap leg, and it is what makes a
  `bool`-returning design visibly worse than an index-returning one.
- ⚠ a field with an entry of exactly `0.0` and one of `-0.0` must be ACCEPTED
  (`-0.0 >= 0.0` is `True` in IEEE, and a naive `np.signbit` implementation would
  reject it). One line, and it is the only place the predicate can be subtly wrong
  without any physics being involved.

### 5.5 What step 4 must NOT be fused with

Memo D §5 measured that the campaign's C1 clause ("DD documented non-cone-preserving
… with step/SC `True`") has **no witness**: `[M]` the scheme registry holds exactly
two schemes, both `is_positivity_preserving = False`; the only `True` values in the
tree are on test doubles. That is a `plan-authoring` §6c gate-without-witness defect
**for the FLAG**, not for the predicate — §5.3 shows the predicate has a real witness
today. ⟹ land the predicate in step 4; leave the flag's `True` arm to whichever step
lands a Step / step-characteristic scheme.

---

## 6. Step-boundary audit (`plan-authoring` §6b/§6c)

### 6.1 ⛔ Step 1 as written creates a transient Pattern-2 TWIN, and 3 committed tests would pin the DEAD copy

The brief's step 1 relocates ρ / true-error / `where_largest` onto the iteration
layer while "the Displacement TYPE still exists after this step". If the three
METHODS also survive on that type, the tree carries two implementations of each and
`tests/numerics/test_affine_flux_algebra.py`'s three diagnostic rows (§3.1) assert
the copy production no longer calls.

`[M]` the cost of not fusing is zero-risk to fix, because the methods have almost no
consumers: `grep -rn "true_error_estimate\|where_largest" orpheus/` outside the
defining module returns **prose only — 0 production call sites**; `contraction_ratio`
has exactly **one** (`iteration.py:792`).

⟹ **Step 1 must delete `Displacement.contraction_ratio`, `.true_error_estimate` and
`.where_largest` in the same commit that lands their iteration-layer replacements**,
and re-point those three tests in that commit. The `Displacement` TYPE survives step
1 as a bare marker + registry; only the diagnostics move. Not a step split — a
one-line widening of step 1's scope.

### 6.2 ⛔ Step 1 retires the leaf FINDER, whose replacement must keep the coupled walk

`_flux_displacement_leaf` (`iteration.py:421-450`) does two jobs: WALK the composite
(`.interior`, else `systems[0]`) and DECIDE it has diagnostics (`hasattr(...,
"contraction_ratio")`). Step 1 destroys the second, and the natural repair — delete
the finder — destroys the first. `[M]` `test_psi_half_coupling.py::test_g_d1_7_…`
exists precisely because the coupled walk once went silent.

⟹ step 1's replacement is a STRUCTURAL walk with no `hasattr` test, and the coupled
gate re-points to it **in the same commit**.

### 6.3 ✅ Step 2 does not cut a call-site set

The type distinction has exactly two consumers — `dsa.py:635` (isinstance) and
`iteration.py:780` (the corrector add) — and both are inside step 2. After step 2
nothing mints a `Displacement` (`FluxRole.__sub__` gone, DSA re-typed), so step 3 is
a pure delete of dead code. Coherent.

### 6.4 ✅ Step 4 has its witness (no §6c violation)

F7: a production-tier violating field exists TODAY, reachable through the public
entry. The gate lands with the case it catches.

### 6.5 ⚠ Step 2 must carry the DSA negative test, not step 4

§4.5's guard test is about the arm step 2 rewrites. Writing it in step 4 leaves an
interval in which the rewritten guard is ungated — the same defect as §6c with a
*guard* in place of a gate.

---

## 7. The carve's own mutation battery (what proves the SUITE, not one gate)

Run scoped, never the full suite. Total ≈ 3 min per arm.

| # | mutation (in-process; NEVER `git checkout`) | must RED |
|---|---|---|
| B1 | `Field.__add__` body → `self.values - other.values` | §4.1 a/c/d/e; the bit-identity wall (all 4 cases) |
| B2 | drop the mesh arm in `_bases.py._check_partner` | §4.3 leg (a), every leaf |
| B3 | `DSACorrection.apply` — delete the interior-type guard | §4.5 leg (a) |
| B4 | `DSACorrection.apply` — return the bulk correction with a ZEROED trace arm | the DSA value case (§1.3) on the REFLECTIVE fixture; **must stay GREEN on the vacuum one** — that asymmetry is the activation proof |
| B5 | the relocated ρ uses the whole composite's norm | §2's pin (`[M]` already measured: 2 failed) |
| B6 | the relocated finder drops the `systems[0]` recursion | `test_psi_half_coupling.py::test_g_d1_7_…` only |
| B7 | **positive control** — `StreamingOperator.apply` returns its input | a large fraction of `tests/sn/`; if it does not, the harness is broken (`vv` #17) |

⚠ Budget the battery off the MUTATED cost, not the green one: a garbaged operator
destroys convergence and a run sized on the 1.6 s baseline will time out.

---

## 8. Refuted / unworkable items in the brief (first-class output)

1. ⛔ **"SI computes the diagnostic trajectory from flat norms" (step 1) is not
   value-neutral as stated.** `[M]` the composite flat norm moves ρ by **4.71e-3**;
   only the INTERIOR block's norm reproduces. The interior-leaf metric-vs-flat swap
   IS neutral (`[M]` ≤1 ULP) — but only because the space is Euclidean today (F8).
2. ⛔ **"exercised by DD's existing negative-flux witness (`TestPositivityFailure`)"
   does not work as written.** `[M]` that test lives at
   `tests/sn/sweep/core/test_diamond.py:793-855` (not "the diamond scheme's test
   file" under `transport/spatial/`), and it asserts on
   `strat.update(...).outgoing_spatial_flux` — a bare ndarray from a single CELL
   VISIT. There is no field, so an element predicate on a `Field` cannot be
   exercised by it without wrapping. F7 supplies a strictly better witness: a
   converged production solve with negative entries, one parameter away from a
   positive sibling.
3. ⛔ **"the ~16 affine raise-sites" and "the 7 files pinning the flux+flux
   TypeError" are two OVERLAPPING counts, not two disjoint work items — and both
   under-report.** Measured with the predicates written out:
   `grep -rn "raises(TypeError" tests/ -A2 | grep "affine\|origin"` → **16 sites in
   10 files**; `grep -rln "affine_combination\|flux + flux\|flux_add_flux" tests/` →
   **8 files** (memo D's narrower predicate found 7). `[M]` the **union is 12 files**,
   of which 1 (`test_ordinate_scan.py`) is a false positive ⟹ **11 real files**.
   4 are raise-site-only (`test_composite.py`, `test_radial_characteristic_field.py`,
   `test_scalar_boundary_flux.py`, `test_radial_characteristic_split_leaves.py`) and
   2 are flux+flux-only (`test_2d_l2_matvec_correctness.py`, `test_ordinate_scan.py`).
   Plan against the union; a plan built on either count alone misses 2–4 files.
4. ⛔ **`tests/sn/sweep/core/test_ordinate_scan.py` is a false positive** — its
   "affine" is the DD recurrence, not the field algebra. Deleting it from the
   migration list is part of the plan.
5. ⚠ **"the torsor-form linearity gates … should re-spell trivially" is right about
   `test_declared_law_is_linear.py` and WRONG about
   `test_2d_l2_matvec_correctness.py`**, which does not merely re-spell: its whole
   `.values`-level workaround DELETES and the gate gets STRONGER (§3.2). Budget it as
   a rewrite, and claim the strengthening.
6. ⚠ **`grep -r torsor` returning only the retirement note is a necessary but weak
   done-when.** It cannot see the four labelled equations, the `affine-*` label
   names, or the prose that argues the doctrine without the word. Add: the four
   labels are retired-or-restated, `grep -rn "affine_combination\|_DISPLACEMENT_CLS\|
   _check_torsor_partner" orpheus/ tests/ docs/` is empty, and
   `nexus dead_references` is clean after the doc pass.

---

## 9. Residual gaps, and the open questions for the user

### 9.1 Gaps this plan leaves OPEN (measured absent, not designed away)

| # | gap | why it is left open |
|---|---|---|
| G1 | The DSA (§1.3) and adjoint (§1.4) baseline `.npy` files are **not captured** — I did not write to `tests/sn/_data/`. | Capture belongs in the carve's first commit at the pre-carve HEAD, together with the case rows, so the `.npy` and the case land atomically. Doing it from a planning dispatch would leave data files with no consumer. |
| G2 | No gate anywhere pins that a **coupled/curvilinear** solve is value-neutral through the RC displacement leaves. The DD snapshots cover `sphere_*` / `cyl_*` at `SAFETY × conv_tol`, and 2 of the `cyl_*` rows already pre-drift. | A dedicated stored-value case for a carrying sphere would close it; it is the only consumer class with no bit-identity wall. |
| G3 | `is_positivity_preserving` still has **0 production readers** (`[M]` memo D §3.1, re-checked at `000cf144`), and `scheme.py:528`'s claim that it "gates negative-flux diagnostics" is aspirational. | The cone predicate (§5) is the natural first reader, but wiring it is a behaviour change, not CS3's neutrality claim. File it. |
| G4 | 4 test modules in the blast radius carry **no V&V marker at all** (§3.5). | Free to fix while editing; not a CS3 correctness item. |
| G5 | `orpheus/sn/solver.py:2960` + issue **#353** rest on a doctrine claim that is already false today and doubly false after CS3. | Belongs to the doc/comment sweep (step 5), but the ISSUE needs re-scoping too — its premise, not just its prose. |

### 9.2 Questions the carve should not start without

- **Q1 (blocking, §2.5).** Is the relocated ρ defined on the **space** norm or the
  **Euclidean** norm? `[M]` they differ by ≤1 ULP today and by **1.12e-3** once CS2
  installs a physical metric. Under "space norm" the §2 pin is a CS3 instrument that
  CS2 legitimately reds and re-derives; under "Euclidean" it is permanent. Both are
  defensible; undecided is not.
- **Q2.** Does the relocated diagnostic surface live on `SourceIteration` (attributes,
  as today) or on `IterationRecord` as a named `StoppingCriterion`-style trajectory?
  §2's `_diagnostics()` is a single function precisely so either answer costs one
  edit — but the answer decides whether `test_flux_displacement_diagnostics.py` and
  `test_psi_half_coupling.py::test_g_d1_7_…` re-point to an object they already hold
  or to one they must now thread out of the record tree.
- **Q3.** Does the cone predicate return a `bool` or the offending INDICES?
  `vv` anti-#14 argues for the structure (it makes its own correctness assertable and
  turns §5.4's sharpest leg from unwritable to one line). Naming it also fixes the
  vocabulary for CS4 and Campaign 2.
- **Q4.** Is `affine_combination`'s content (the relaxation blend) kept as a named
  operation anywhere, or does it dissolve into `0.7·ψ₂ + 0.3·ψ₁`? `[M]` 0 production
  callers, 5 test consumers — dissolving is defensible, but the RELAXATION concept is
  real and a future consumer will want a name for it (`coding-elegance` Pattern 3).

### 9.3 Run recipes (scoped; never the full suite)

```
# the new capture gate                         [M] 5 passed, 0.46 s
#   (whole dir clean too: tests/numerics/ [M] 2349 passed, 5:32)
.venv/bin/python -O -m pytest tests/numerics/test_si_diagnostic_trajectory.py -q

# the bit-identity wall (escalated)            [M] 3 passed, 1.60 s
.venv/bin/python -O -m pytest tests/sn/solve/test_affine_carve_bit_identity.py -q \
  -W error::tests.sn.regression._regression_assert.DriftWarning

# the broad value wall (2 characterised pre-drifts deselected)   [M] 11 passed, 52.7 s
.venv/bin/python -O -m pytest tests/sn/regression/test_dd_regression.py -q \
  -W error::tests.sn.regression._regression_assert.DriftWarning \
  --deselect "tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_1g_homogeneous_folded_2x4_dd_n20]" \
  --deselect "tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_1g_homogeneous_folded_4x8_dd_n20]"

# the algebra + diagnostics migration surface (the 11-file union of §8.3)
.venv/bin/python -O -m pytest tests/numerics/test_affine_flux_algebra.py \
  tests/transport/ tests/sn/primitives/ tests/sn/acceleration/ \
  tests/sn/solve/test_flux_displacement_diagnostics.py -q
```
