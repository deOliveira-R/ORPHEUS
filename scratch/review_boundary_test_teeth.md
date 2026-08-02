# Boundary test suite — does it have TEETH?

QA quadrant of the boundary-machinery complete-scrutiny review.
Branch `refactor/operator-strategy-layers` @ `73627b71`. Host `.venv` (Py 3.14),
canonical invocation `.venv/bin/python -O -m pytest`, SERIAL.

Evidence key: **[M]** measured this session · **[R]** read (quoted + `file:line`)
· **[G]** grep/exhaustive · **[U]** UNVERIFIED (may not enter the plan).

---

## HEADLINE — the Mode-8 hypothesis is REFUTED

> **The fraction of the boundary suite's assertions that are inert under the
> canonical `python -O` invocation is `0 / 676 = 0.0 %`.**

The brief's premise ("`-O` strips bare `assert` to no-ops, therefore the
bare-assert-heavy `tests/geometry/` files are inert") is **false for this
codebase**, because pytest's assertion **rewriter** replaces `assert` statements
with explicit `if not …: raise AssertionError(…)` in every *collected* test
module — before the interpreter ever sees them. `-O` only strips asserts in
modules the rewriter does **not** touch (production `orpheus/`, non-collected
helper scripts).

61.7 % of the suite's assertions ARE bare `assert` — that number is real and is
reported below — but **all 417 of them live in collected `tests/` modules and
all 417 fire under `-O`.** Measured two ways, below.

**This is not a clean bill of health.** The suite's real weakness is not Mode 8;
it is Mode 11/12 (§4–§7): a large majority of the bare asserts are
`isinstance`/`is`/`==`-on-a-tag *structural smoke*, not value gates — and the
seven laws are gated with radically uneven force.

---

## 1. Mode-8 `-O` exposure — MEASURED

### 1.1 Synthetic control [M]

`$CLAUDE_JOB_DIR/tmp/test_mode8_probe.py`, three deliberately-failing tests
(bare `assert`, `np.testing.assert_allclose`, `pytest.fail`):

| invocation | result |
|---|---|
| `.venv/bin/python -m pytest` (no `-O`) | **3 failed** in 0.16 s |
| `.venv/bin/python -O -m pytest` | **3 failed** in 0.09 s (+1 warning) |

The bare-assert test failed **identically** under `-O`. The extra warning is
pytest's own `PytestConfigWarning: assertions not in test modules or plugins
will be ignored` — note the wording: *not in test modules*. It is a warning
about production/helper modules, **not** about the collected tests.

### 1.2 Realistic probe on two ACTUAL boundary test files [M]

Copies (never the repo files) placed at
`$CLAUDE_JOB_DIR/tmp/realprobe/`, one bare assert falsified in each:

* `test_bc_universal_invariants.py:399` `assert isinstance(bc.source, NoSource)`
  → `… and False, 'MUT-A bare assert must FAIL'`
* `test_law_composition.py:69` `assert isinstance(tree, LawSum)`
  → `… and False, 'MUT-B bare assert must FAIL'`

| invocation | result |
|---|---|
| no `-O` | `2 failed, 63 passed` — `MUT-A`, `MUT-B` both RED |
| **`-O` (canonical)** | `2 failed, 63 passed` — `MUT-A`, `MUT-B` both RED |

**Byte-identical verdict in both modes.** Bare asserts in `tests/` are teeth.

### 1.3 Where Mode-8 exposure COULD live — exhaustive [G]

| location | bare `assert` count | `-O` status |
|---|---|---|
| `orpheus/geometry/boundary/**` | **0** | n/a |
| `orpheus/sn/boundary/**` | **0** | n/a |
| `orpheus/sn/operators/boundary.py` | **0** | n/a |
| `orpheus/diffusion/**` | **0** | n/a |
| `orpheus/numerics/operator.py` | **0** | n/a |
| non-collected helpers under `tests/` | **1** | INERT — see below |

The single inert bare assert **[R]**
`tests/geometry/_generate_bc_equivalence_snapshots.py:339`:

```python
assert case.case_id == "mixed_30spec_70white_LS4", (
```

**[G]** That module is **never imported by a test** — the only references
(`test_bc_equivalence_snapshot.py:17,49,79`) are docstring prose and an error
message telling a human to run `python -m tests.geometry._generate_bc_equivalence_snapshots`.
It is a manually-invoked baseline generator. Its assert is a developer guard,
not a suite gate, and it is inert under `-O` — a real but *low-severity*
Mode-8 site (severity: it silently permits regenerating a snapshot under the
wrong case-id if a future author runs the generator with `-O`).

### 1.4 The census [M] — AST-classified, all 20 boundary-related files

`python -m ast` over the file set (`$CLAUDE_JOB_DIR/tmp/count_asserts.py`):

| file | test fns | collected | bare | np.testing | raises/warns | pytest.fail | `assert_*` calls |
|---|---:|---:|---:|---:|---:|---:|---:|
| `tests/geometry/test_bc_universal_invariants.py` | 47 | 47 | 16 | 9 | 14 | 0 | 32 |
| `tests/geometry/test_boundary_trace_law.py` | 14 | 14 | 31 | 0 | 1 | 0 | 6 |
| `tests/geometry/test_bc_errors.py` | 11 | 11 | 35 | 0 | 9 | 0 | 0 |
| `tests/geometry/test_boundary.py` | 25 | 25 | 18 | 18 | 3 | 0 | 0 |
| `tests/geometry/test_bc_equivalence_snapshot.py` | 8 | 8 | **0** | 9 | 0 | 0 | 0 |
| `tests/geometry/test_law_composition.py` | 18 | 18 | 56 | 4 | 2 | 0 | 0 |
| `tests/geometry/test_reduced_operator.py` | 31 | 47 | 78 | 0 | 5 | 0 | 0 |
| `tests/sn/operators/test_sn_boundary_realizer.py` | 25 | 25 | 42 | 12 | 6 | 0 | 0 |
| `tests/sn/operators/test_sn_boundary_operator.py` | 11 | 20 | 15 | 6 | 2 | 0 | 0 |
| `tests/sn/operators/test_boundary_conditions.py` | 11 | 11 | 17 | 0 | 2 | 0 | 0 |
| `tests/sn/operators/test_snmesh_realizer_wiring.py` | 11 | 11 | 26 | 12 | 3 | 3 | 0 |
| `tests/sn/operators/test_angular_average_operator.py` | 16 | 16 | 8 | 8 | 4 | 0 | 0 |
| `tests/sn/operators/test_bc_extraction_2d.py` | 6 | 8 | 7 | 4 | 0 | 0 | 0 |
| `tests/sn/operators/test_bc_extraction_matvec.py` | 7 | 33 | 5 | 4 | 0 | 0 | 1 |
| `tests/diffusion/test_boundary_realizer.py` | 21 | 35 | 12 | 7 | 6 | 0 | 0 |
| `tests/numerics/test_periodic_wrap_operator.py` | 5 | 5 | 5 | 6 | 0 | 0 | 0 |
| `tests/transport/fields/test_angular_boundary_flux.py` | 36 | 36 | 25 | 17 | 12 | 0 | 0 |
| `tests/transport/fields/test_scalar_boundary_flux.py` | 15 | 15 | **0** | 14 | 8 | 10 | 0 |
| `tests/sn/primitives/test_boundary_face_layout.py` | 5 | 5 | 20 | 0 | 0 | 0 | 0 |
| `tests/geometry/_generate_bc_equivalence_snapshots.py` | 0 | **not collected** | 1 | 0 | 0 | 0 | 0 |
| **TOTAL** | **323** | **390** | **417** | **130** | **77** | **13** | **39** |

* total assertion-bearing statements **676**
* bare-assert fraction **417/676 = 61.7 %**
* **inert-under-`-O` fraction 0/676 = 0.0 %** (the 1 non-collected assert is
  outside the suite entirely; including it: 1/677 = 0.15 %)

`assert_*` calls (39) are the harness/law invariant methods
(`assert_is_involutive`, `assert_submarkov`, `assert_realizable`, …) — these are
**function calls**, immune to `-O` regardless of whether their *bodies* use bare
asserts (they do not: 0 bare asserts in `orpheus/geometry/boundary/**`; they
`raise` typed errors).

**Conclusion for §8.4 of the review plan:** delete the Mode-8 concern. The
correct concern for this subsystem is what the 417 bare asserts *assert*
(§4), not whether they run.

---

## 2. The boundary test surface — enumerated

**[G]** Files that **import** `orpheus.geometry.boundary` or
`orpheus.sn.boundary` (the machinery-under-test set), plus files that consume a
realized boundary operator by name:

| tier | file | collected |
|---|---|---:|
| **law layer** | `tests/geometry/test_boundary_trace_law.py` | 14 |
| | `tests/geometry/test_bc_errors.py` | 11 |
| | `tests/geometry/test_bc_universal_invariants.py` | 47 |
| | `tests/geometry/test_law_composition.py` | 18 |
| | `tests/geometry/test_bound_compat.py` | *(in set; see below)* |
| **realizer (SN)** | `tests/sn/operators/test_sn_boundary_realizer.py` | 25 |
| | `tests/sn/operators/test_snmesh_realizer_wiring.py` | 11 |
| | `tests/sn/operators/test_boundary_conditions.py` | 11 |
| | `tests/sn/operators/test_angular_average_operator.py` | 16 |
| | `tests/geometry/test_boundary.py` | 25 |
| **realizer (diffusion)** | `tests/diffusion/test_boundary_realizer.py` | 35 |
| **snapshot** | `tests/geometry/test_bc_equivalence_snapshot.py` | 8 |
| **composite operator** | `tests/sn/operators/test_sn_boundary_operator.py` | 20 |
| | `tests/sn/operators/test_bc_extraction_2d.py` | 8 |
| | `tests/sn/operators/test_bc_extraction_matvec.py` | 33 |
| **numerics primitives** | `tests/numerics/test_periodic_wrap_operator.py` | 5 |
| | `tests/numerics/test_incoming_ordinate_mask_tensor.py` | *(in set)* |
| | `tests/geometry/test_reduced_operator.py` | 47 |
| **field carriers** | `tests/transport/fields/test_angular_boundary_flux.py` | 36 |
| | `tests/transport/fields/test_scalar_boundary_flux.py` | 15 |
| | `tests/sn/primitives/test_boundary_face_layout.py` | 5 |
| **generator (NOT collected)** | `tests/geometry/_generate_bc_equivalence_snapshots.py` | 0 |

**[M]** The 18-file mutation-harness set collects **323 passed, 3 skipped** at
baseline under `-O` in ≈2.3 s.

Additional consumers found **[G]** but not in the harness set (they exercise a
realized boundary operator inside a bigger solve):
`tests/sn/architecture/test_monomorphic_leaves.py`,
`test_stage_separation.py`, `tests/sn/operators/test_capability_survival.py`,
`test_g_adjoint_reciprocity.py`, `test_inverse_adjoint_coherence.py`,
`test_ld_adjoint_deferral.py`, `test_operator_block_role.py`,
`test_apply_full_field_codomain.py`, `test_one_representation_instance.py`,
`tests/sn/solve/test_d3_admission.py`, `test_gauss_seidel_reification.py`,
`tests/sn/sweep/core/test_phase_c_gates.py`,
`tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py`,
`tests/sn/eigenvalue/test_keff_2d.py`,
`tests/sn/verification/analytical/test_prescribed_inflow_consistency.py`,
`tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py`,
`tests/transport/test_method.py`, `tests/numerics/test_inverse_universal.py`,
`tests/numerics/test_registry_mixin.py`,
`tests/numerics/test_operator_capability_predicates.py`,
`tests/sn/primitives/test_method_space.py`,
`tests/diffusion/test_operators.py`.

---

## 3. Mutation sweep — WHICH claims have teeth [M]

All mutations are **in-process monkeypatches** installed by a throwaway pytest
plugin (`$CLAUDE_JOB_DIR/tmp/mutplugin.py`, `mutplugin2.py`); **no file under
`orpheus/` or `tests/` was written**, and every mutation ships a **bite check**
that aborts the run if the production behaviour did not actually change.

> **Self-catch, recorded because it is the exact trap this review hunts.**
> My FIRST albedo mutation installed `AlbedoBoundary.__post_init__`. That is a
> **no-op**: a dataclass bakes the `self.__post_init__()` call into `__init__`
> only if the hook existed at *decoration* time, and **[M]**
> `hasattr(AlbedoBoundary, "__post_init__") == False`. The run reported
> `303 passed` — i.e. **"no test catches an albedo inversion"** — which was a
> **false finding produced by a vacuous mutation**. Re-done at the realizer
> seam with a bite check, the same mutation reddens **7** tests. Every
> mutation below carries a bite check for this reason.

### 3.1 Leaf-action mutations (the realized operator's numbers)

| # | mutation | failure mode | tests RED | verdict |
|---|---|---|---:|---|
| M1 | reflective: permutation → **identity** | Mode 2 | **13** | GATED |
| M2 | reflective: permutation → `np.roll(perm,1)` | Mode 5 | **12** | GATED |
| M3 | white: drop the `\|Ω·n\|` cosine factor | Mode 3 | **11** | GATED |
| M4 | white: normalise by `Σw` not `Σ w\|μ\|` | Mode 3 | **16** | GATED |
| M5 | white: no outgoing-hemisphere mask | Mode 6 | **16** | GATED |
| M6 | vacuum: mask the **outflow** rows instead | Mode 6 (ERR-041 shape) | **8** | GATED |
| M7 | prescribed_inflow: drop the `Γ₋` mask | ERR-047 re-drop | **3** | GATED |
| M8 | `_reflect_trace`: drop the `sel` row restriction | — | **5** | GATED |
| M9 | SN α → 1−α (albedo/reflective/white) | Mode 1/3 | **7** | GATED |
| M10 | diffusion α → 1−α | Mode 1/3 | **3** | GATED |
| M11 | diffusion `zero_flux`: 𝒜 = −1 → **+1** | Mode 1 (sign flip) | **4** | GATED |
| M12 | diffusion `vacuum`: 𝒜 = 0 → **1** | Mode 6 | **3** | GATED |

Named catchers worth recording:

* M1/M2 → `test_sn_boundary_realizer.py::TestRealizeReflective::test_specular_unit_albedo_lebedev_matches_hand_computed`
  and `…::test_specular_partial_albedo_matches_hand_computed` (**hand-computed**
  references — genuine L0 value gates, not smoke).
* M3/M4/M5 → `test_angular_average_operator.py::TestCosineWeightedCurrentConservation::*`
  **and** `TestSelfAdjointnessCosineWeighted::*` — the white conservation claim
  is a REAL gate (§5).
* M6 → `test_sn_boundary_realizer.py::TestRealizeVacuum::test_lebedev17_{xmin,ymax}_zeroes_only_inflow`.
* M7 → `test_bc_universal_invariants.py::TestSourceLivesOnIncomingTraceInvariant::test_delivered_q_vanishes_off_the_incoming_trace`.
* M8 → `test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[…]`
  ×4 geometries, failing with *"B emitted non-zero on the outflow rows"*.

> **[M] Cross-reference to review-plan §4.1.** The universal discard of the
> outflow rows in `_reflect_trace` is **not** an accident awaiting a consumer —
> it is an **asserted contract**: M8 makes `B` preserve the outflow rows and
> **four** tests go red *demanding* they be zero. Any plan step that
> "restores outflow preservation" must re-baseline those four assertions and
> the 2-D balance gate. That is a design decision with a committed test, not a
> latent gap.

> **[M] A Mode-12 observation worth the plan's attention.** The physics gate
> `test_bc_extraction_2d.py::…::test_reflective_boundary_balance_at_convergence`
> (2-D reflective ⇒ `keff == k_inf`) catches **M1** (identity "reflection"
> — k moves to 1.5629 vs k_inf) but is **BLIND to M2** (permutation rolled by
> one). A rolled permutation is still a permutation, so total in = total out
> and the balance functional annihilates the index drift. The catchers for
> Mode-5 index drift are the **hand-computed structural gates**, never the
> eigenvalue. Do not let a plan cite the keff gate as reflection-table coverage.

### 3.2 Guard-disabling mutations (does each production raise have a catcher?)

Every law invariant / realizer refusal was neutered one at a time.
**Every single one reddens exactly its own negative test.** No dead guard.

| guard neutered | tests RED | catcher |
|---|---:|---|
| `ReflectiveBoundary.assert_is_involutive` (ERR-044) | 1 | `TestReflectiveInvolutionInvariant::test_raises_for_non_involutive_perm` |
| `…assert_geometry_map_measure_preserving` (ERR-042) | 2 | `TestReflectiveMeasurePreservingInvariant::{test_raises_for_involutive_weight_class_mispairing, test_raises_at_realize_time}` |
| `…assert_reflection_maps_inflow_to_outflow` (ERR-045) | 2 | `TestReflectiveInflowToOutflowInvariant::*` |
| `WhiteBoundary.assert_submarkov` (ERR-046) | 1 | `TestWhiteSubmarkovInvariant::test_raises_for_supercritical_albedo` |
| `WhiteBoundary.assert_response_positive_if_declared` (ERR-043) | 1 | `TestWhiteResponsePositiveInvariant::test_raises_for_negative_albedo` |
| `AlbedoBoundary.assert_submarkov` (ERR-046) | 1 | `TestAlbedoSubmarkovInvariant::…` |
| `AlbedoBoundary.assert_response_positive_if_declared` (ERR-043) | 1 | `TestAlbedoResponsePositiveInvariant::…` |
| `BoundaryTraceLaw.assert_source_lives_on_incoming_trace` (ERR-047) | 2 | `TestSourceLivesOnIncomingTraceInvariant::*` |
| SN `ZeroFluxBoundary` refusal | 1 | `test_boundary.py::test_sn_realizer_refuses_zero_flux` |
| SN ERR-041 vacuum orientation guard | 2 | `TestVacuumTraceOrientationGuard::{test_swapped_annotation_raises, test_single_contaminated_index_raises}` |
| diffusion `PeriodicBoundary` refusal | 2 | `TestRefusals::test_periodic_refused` + `TestComposition::test_composition_refusal_propagates` |
| diffusion `PrescribedInflow` refusal | 1 | `TestRefusals::test_prescribed_inflow_refused` |
| `realize_recursively`: `LawScaled` scalar dropped | **7** | 5 in `test_law_composition.py` + 2 in `test_boundary_realizer.py` |
| `realize_recursively`: `LawSum` right operand dropped | **5** | 4 in `test_law_composition.py` + 1 diffusion |

**This is the strongest single result of the review: the boundary subsystem's
guards are, without exception, mutation-live.** That is a materially better
posture than the declared-but-unconsumed *capabilities* the review found on the
production side.

---

## 4. The one genuinely toothless file — `tests/geometry/test_bc_errors.py`

**[R]** All 11 tests share this shape (`test_bc_errors.py:47-61`, representative):

```python
err = BoundaryError("msg", law="my_law")
assert isinstance(err, ValueError)
assert isinstance(err, BoundaryError)
assert err.law == "my_law"
assert str(err) == "msg"
with pytest.raises(BoundaryError):
    raise err
```

**[G]** Every one of the file's 9 `pytest.raises` blocks contains exactly
`raise err`, where `err` was constructed two lines above as an instance of the
very class being caught.

**`with pytest.raises(X): raise err` where `isinstance(err, X)` is a
TAUTOLOGY** — the `vv-principles` Mode-8 "fires-but-cannot-fail" class. **No
input exists that makes it red.** It cannot detect: a production site that
stops raising, a production site that raises the wrong subclass, a changed
message, or a removed `law=` attribute at any raise site.

What the file *does* legitimately pin: the class hierarchy
(`BoundaryError <: ValueError`, each typed error `<: BoundaryError`) and the
`law=`/`str()` attribute contract. Those parts are real (`assert isinstance(…)`
would red if the base class changed) — **only the `pytest.raises` legs are
tautological.**

**[M]** Confirmed by measurement: **not one** of the 14 guard-disabling
mutations in §3.2 reddened a single test in `test_bc_errors.py`
(11/11 green throughout). Every production raise site is covered by
`test_bc_universal_invariants.py` / `test_sn_boundary_realizer.py` /
`test_boundary.py` / `test_boundary_realizer.py` instead.

**Verdict:** not a coverage HOLE (the raises are covered elsewhere) but a
**coverage-count inflation**: 11 tests + 35 bare asserts + 9 `pytest.raises`
that read as "the error paths are heavily tested" while contributing zero
mutation-detection. Recommend: delete the 9 tautological
`with pytest.raises(...): raise err` legs, keep the hierarchy/attribute
assertions, or fold the whole file into a single parametrized hierarchy test.


---

## 5. The seven laws × two methods — reachability and gating [M]

### 5.1 Ground truth (realized in-process, `Quadrature.lebedev(5)`, face `xmax`)

| law | SN registry | diffusion registry | SN realized | diffusion realized (𝒜) |
|---|:--:|:--:|---|---|
| `vacuum` | **yes** | **yes** | `TensorProductOperator` (mask ⊗ I) | `ZeroOperator` (𝒜 = 0) |
| `reflective` | **yes** | **yes** | `TensorProductOperator` (perm ⊗ I) | `IdentityOperator` (𝒜 = 1) |
| `white` | NO | NO | `TensorProductOperator` (avg ⊗ I) | `IdentityOperator` (𝒜 = 1) |
| `albedo` | NO | **yes** | `ScaledOperator` | `ScaledOperator` (𝒜 = 0.37) |
| `periodic` | NO | NO | `TensorProductOperator` — **acts as the IDENTITY** | **REFUSES** (`BoundaryError`) |
| `prescribed_inflow` | NO | NO | `IncomingSourceOperator` | **REFUSES** (`BoundaryError`) |
| `zero_flux` | NO | **yes** | **REFUSES** (`BoundaryError`) | `ScaledOperator` (𝒜 = −1) |

**[R]** Confirms review-plan §6 (was `[U]`): `SNMesh.BOUNDARY_OPERATOR_REGISTRY`
= `{"vacuum": VacuumInflow, "reflective": ReflectiveBoundary}`
(`orpheus/sn/mesh/augmented_mesh.py:171-173`);
`DiffusionMesh.BOUNDARY_OPERATOR_REGISTRY` =
`{"vacuum", "reflective", "albedo", "zero_flux"}`
(`orpheus/diffusion/augmented_mesh.py:158-163`).
**5 of 7 laws are unreachable from an SN `BC` tag; 3 of 7 from a diffusion tag.**

**[M] New instance of the review's "declared capability, no consumer" pattern —
in its purest form.** `PeriodicBoundary` realizes for SN to
`PeriodicWrapOperator() & IdentityOperator()`, and **[M]** its apply is
`np.array_equal(op.apply(x), x) == True` — the **identity**. **[R]**
`tests/numerics/test_periodic_wrap_operator.py` (5 tests) gates precisely that
it IS the identity (`test_apply_is_identity_various_shapes`,
`test_apply_transpose_is_identity`, `test_apply_returns_fresh_copy`), and
**[R]** the realizer comment concedes it: *"the current `PeriodicWrapOperator`
body is identity-with-copy (the SN sweep handles the spatial wrap via its
face-pair indexing)… When `PeriodicWrapOperator` gains a non-trivial spatial
pushforward (follow-up issue), the second factor will carry that structure"*
(`orpheus/sn/boundary/realizer.py:324-331`). So periodic is: **not declarable,
realizes to a no-op, and is fully tested — for being a no-op.** Its 5 tests are
honest about what they gate; the *capability* is the empty one.

### 5.2 GATED vs SMOKE, per law × method [M — from the §3 mutation sweep]

`V` = a mutation of the realized ACTION reddens a test.
`E` = error-path (refusal / invariant) gated, mutation-verified.
`S` = structure/type/shape smoke only.

| law | SN | diffusion |
|---|---|---|
| `vacuum` | **V** (M6, 8 red) + **E** (ERR-041 orientation guard, 2 red) | **V** (M12, 3 red) |
| `reflective` | **V** (M1 13 red / M2 12 red / M9 α) + **E** (3 table invariants ERR-042/044/045, 5 red) | **V** (M10, 3 red) |
| `white` | **V** (M3 11 / M4 16 / M5 16 / M9 α) + **E** (ERR-043/046, 2 red) | **V** (M10, shares the α leg) |
| `albedo` | **V** (M9, 7 red) + **E** (ERR-043/046, 2 red) | **V** (M10, 3 red) |
| `periodic` | **S only** — the action IS the identity; nothing numeric to gate | **E** (refusal, 2 red) |
| `prescribed_inflow` | **V** (M7, 3 red — the ERR-047 Γ₋ mask) + **E** (2 red) | **E** (refusal, 1 red) |
| `zero_flux` | **E** (refusal, 1 red) | **V** (M11 sign flip, 4 red) |

**Every law × method cell is either action-gated or refusal-gated, and every
one was mutation-verified.** The only `S`-only cell is SN-periodic, and it is
`S`-only *because there is no action*.

---

## 6. Is the white BC's conservation claim tested? YES — and it has teeth [M]

**[R]** The claim, `orpheus/sn/boundary/angular.py:55-58`:

> "The denominator is the outgoing cosine-weighted weight sum — this
> normalization makes the operator conservative for any quadrature: the
> cosine-weighted outgoing current equals the cosine-weighted incoming current
> (per Bell & Glasstone 1970 §1.5)."

**[R]** The gate:
`tests/sn/operators/test_angular_average_operator.py:56-94`,
`@pytest.mark.l1 class TestCosineWeightedCurrentConservation` — two tests
(`test_outgoing_current_equals_incoming_lebedev` on `lebedev(17)` with a random
positive ψ; `test_outgoing_current_with_anisotropic_input` on
`level_symmetric(6)` with `ψ = exp(μ_z)(μ_z + 2)`), each asserting
`j_in == j_out` at `rtol=1e-13` via `np.testing.assert_allclose`.

**[M] Mutation-verified — RED under all three white mutations** (M3 drop-cosine,
M4 wrong-norm, M5 no-hemisphere-mask). Companion
`TestSelfAdjointnessCosineWeighted` (⟨Ax,y⟩_w = ⟨x,Ay⟩_w) reds under M3 and M5.

**Honest scope, three caveats the plan should carry:**

1. **Procedural, not structural, independence (L11).** The test re-encodes the
   production weight formula: `cos_w_out = quad.weights * quad.mu_x` masked to
   `mu_x > 0` — literally `angular.py:181-184`. Algebraically the identity
   `⟨cos_w, Aψ⟩ = ⟨cos_w, ψ⟩` reduces to `norm == cos_w.sum()`, so the test
   *directly* pins the normalisation and pins the weight formula only
   *indirectly* (via the test-side/production-side mismatch a mutation
   creates). That indirect pin is real — measured — but it is a
   re-transcription, not an independent derivation. A structurally-independent
   reference would be a closed-form `∫_{hemisphere} |Ω·n̂| dΩ = π` check, which
   the suite does **not** have. **[G]**
2. **Operator-level only.** White is not in `SNMesh.BOUNDARY_OPERATOR_REGISTRY`
   (§5.1), so **[G]** no end-to-end SN solve exercises this conservation
   property. The claim is verified for the *primitive*, never for the *solver*.
3. **The two tests the docstrings call "hand-computed" are NOT catchers — [M].**
   `test_boundary.py::test_white_bc_4_point_quadrature_hand_computed`
   (docstring: *"White BC on 4-point GL: explicit hand calculation"*) and
   `::test_white_bc_axis_z_on_product_quadrature` feed a
   **hemisphere-CONSTANT** ψ (1.0 outgoing / 7.0 incoming; 2.0 / 9.0) and assert
   the output equals that constant. A normalised weighted average of a constant
   is that constant **for any weights**, so these tests measure a functional
   whose invariance group contains the entire weight formula (`vv-principles`
   Mode 12). Measured, in isolation:

   | test | control | under M3 (`|μ|` dropped) |
   |---|---|---|
   | `test_white_bc_4_point_quadrature_hand_computed` | PASSED | **PASSED** ← blind |
   | `test_white_bc_axis_z_on_product_quadrature` | PASSED | **PASSED** ← blind |
   | `test_white_bc_returns_cosine_weighted_average` | PASSED | **FAILED** ← the catcher |

   There is no hand-computed number anywhere in the first test — the expected
   `1.0` falls out of normalisation alone. **The docstring over-claims.** The
   real catchers are `test_white_bc_returns_cosine_weighted_average`, the two
   conservation tests, the two self-adjointness tests, and the three white
   snapshots.

---

## 7. The rank-N composition path — `tests/geometry/test_law_composition.py`

**[R] It gates BOTH structure and values — but the value reference is not
independent of the leaves.**

* **Structure (12 tests, `@pytest.mark.foundation`)**: algebra closure
  (`LawSum`/`LawScaled`, constant folding, `/`, unary `-`, non-flattening),
  the "descriptors carry no `apply`" invariant, walker type-family preservation,
  `TypeError` on a non-law node. All `isinstance`/`is`/`==` — structural smoke
  by construction, and appropriately so.
* **Values (2 tests, `@pytest.mark.l1`)**:
  `test_realize_recursively_apply_matches_pointwise_weighted_sum` (the Marshak
  `0.3·Reflective + 0.7·White` case, `lebedev(17)`, `ψ` shape `(N,5,3)`) and
  `test_realize_recursively_nested_apply_matches_distributive_form`
  (`0.5·(0.3·spec + 0.7·white)`), both `assert_array_almost_equal_nulp(nulp=4)`.

**Answer to the review-plan question "is the Marshak case checked numerically
against anything independent?" — NO, and by design.** **[R]** `:277`:

```python
expected = 0.3 * spec_realised.apply(psi) + 0.7 * white_realised.apply(psi)
```

`spec_realised` / `white_realised` are the **production leaves**. The reference
is independent of the **walker**, not of the **leaves**. The test's own
docstring is accurate about this (*"the explicit pointwise weighted sum … uses
only numpy addition and scalar multiplication"*) — the independence claimed is
for the *composition algebra*, and that claim is true.

The subsystem's composite correctness therefore decomposes as
**(leaf value gates, §3.1) ∘ (walker distributivity gate, here)**, which is a
legitimate factorisation — provided the plan does not cite the composition
tests as evidence for leaf values. If both leaves were wrong identically, these
two tests stay green.

**[M] Both value tests are mutation-live:**

| walker mutation | tests RED |
|---|---:|
| `LawScaled` scalar dropped | **7** (5 in `test_law_composition.py`, 2 in `tests/diffusion/test_boundary_realizer.py::TestComposition`) |
| `LawSum` right operand dropped | **5** (4 + 1) |

**[G] Production reach of rank-N**: `realize_recursively` has **no production
caller that builds a `LawSum`/`LawScaled`** — no `BOUNDARY_OPERATOR_REGISTRY`
tag maps to a composite, and the Marshak `0.3·spec + 0.7·white` form appears
only in tests. Rank-N is a *tested-but-unreachable* capability: a fourth
instance of the review's declared-capability pattern, but unlike the other
three it is **fully mutation-gated**, so materialising a consumer would land on
verified machinery.

---

## 8. Adjoint coverage [M]

**[R] The "refusal" is a TWO-part contract, not one.** `is_adjointable` is a
**declared predicate** defaulting to `False` on `LinearOperator`
(`orpheus/numerics/operator.py:623-634`: *"Default `False`; an operator with a
working `apply_transpose` overrides"*), and `adjointable(op)` is just
`return op.is_adjointable` (`:1113-1125`). So `AngularAverageOperator` refuses
by (a) not defining `apply_transpose` **and** (b) not overriding
`is_adjointable`. *My first mutation added only (a) and left the refusal intact
— worth stating because a plan step that "adds `apply_transpose` to white"
without lifting the predicate changes nothing.*

| # | adjoint mutation | tests RED | catcher |
|---|---|---:|---|
| A1 | `SNBoundaryOperator.is_adjointable`: per-face `all` → `any` | **1** | `test_sn_boundary_operator.py::TestApplyTransposeCapability::test_adjointability_drops_when_a_face_lacks_it` |
| A2 | `AngularAverageOperator` advertises the unweighted transpose (both halves) | **1** | `test_angular_average_operator.py::TestPredicates::test_apply_only` |
| A3 | `_reflect_trace` transpose: drop the codomain-row mask | **9** | 4× `TestApplyTransposeCapability::test_adjointable_when_all_faces_support[…]` + 4× `test_monomorphic_leaves.py::test_hilbert_adjoint_reciprocity[B-{slab,sphere,cylinder,cart2d}]` + `test_a_globally_constant_metric_makes_reciprocity_blind` |
| A4 | ray-corner transpose mirror swapped (`RadialCharacteristicBoundaryOperator`) | **3** | `test_psi_half_coupling.py::TestB_b_RayBoundary::{test_transpose_is_exact_euclidean_mirror, test_euclidean_transpose_is_the_vcell_hilbert_adjoint}` + `TestCoupledBuilder::test_forward_block_adjoint_reciprocity` |
| A5 | `_RULED_CORNER_KINDS` widened to admit white/albedo/periodic | **1** | `test_psi_half_coupling.py::TestB_b_RayBoundary::test_unruled_outer_law_is_loud_deferred` |

**Answers to the brief's §7:**

* **Is the adjoint path gated?** YES, and strongly — A3's 9 reds include the
  G-metric Hilbert-adjoint reciprocity gate on all four geometries, with its
  `test_a_globally_constant_metric_makes_reciprocity_blind` Mode-12 control leg
  live and reddening alongside.
* **Is the REFUSAL itself gated?** YES for both refusals, but by **exactly one
  test each** (A2 for white's `apply_only`; A5 for the ray-corner
  `NotImplementedError`). Thin, but live and mutation-verified.
* **[G] One ungated defensive raise.** `SNBoundaryOperator._reflect_trace`
  raises `MissingAdjoint` per face when a caller bypasses `is_adjointable`
  (`orpheus/sn/operators/boundary.py:~310`). Grep across all of `tests/` finds
  `MissingAdjoint` only in `tests/numerics/test_coupled_operator.py` and
  `tests/numerics/test_outer_dyad.py` — **no boundary test drives it**. The
  production comment concedes the status: *"unreachable-in-practice because
  `is_adjointable` gates the composite eagerly, but the per-face raise keeps the
  refusal loud if a caller bypasses the predicate."* Accurate self-description;
  low severity; note it rather than test it.

---

## 9. The toothless / at-risk inventory — six items, all measured

Ordered by how badly the artefact's *label* over-states its *content*.

### 9.1 `catches("ERR-052")` is a coverage claim with no teeth **[M]**

**[G]** `tests/sn/operators/test_boundary_conditions.py:145` carries
`@pytest.mark.catches("ERR-052")` and is the **only** test in the repo with
that marker. **[R]** Its docstring: *"Catches ERR-052: power-iteration flux
non-renormalisation. Without the per-step `ψ /= ‖ψ‖` projection … subcritical
cases … underflow to denormalised FP within ~30-60 outer iterations."*

**[M]** I re-introduced the EXACT documented bug in-process by making the
`isinstance(solver, ProductionRateSolver)` gate at
`orpheus/numerics/eigenvalue.py:255` never match — the pre-fix un-normalised
trajectory, verbatim. **Bite-verified**, not assumed:

| | `renorm_calls` | `\|φ\|max` (vacuum) | `\|φ\|min` | `n_outer` | `k_vac` |
|---|---:|---:|---:|---:|---|
| baseline | **6** | 7.5955e+00 | 1.0381e+00 | 6 | 0.179044149673 |
| bug re-introduced | **0** | 6.1197e-01 | 8.3640e-02 | 6 | 0.179044144225 |

The renormalisation is genuinely gone (12× flux-magnitude shift), and the test
**still passes**: `pytest tests/sn/operators/test_boundary_conditions.py` gives
`11 passed` with and without the mutation.

**Why it went blind:** the case converges in **6** outer iterations; ERR-052's
denormalisation needs **30–60**. And the assertion is an ORDERING
(`k_refl > k_vac`, i.e. `1.875 > 0.179` — a 10× margin), which survives even a
badly degraded iterate. **The marker was a true catcher when written and the
configuration has since drifted out of the failure regime.** Per
`vv-principles` §"a `catches` marker is a COVERAGE CLAIM", this marker is now a
phantom edge in the ERR coverage audit. Fix: either drive the test into the
30–60-outer regime (tighter outer tolerance / a genuinely subcritical
multiplying config) and re-verify the red, or move the marker to a test that
does redden.

### 9.2 `tests/geometry/test_bc_errors.py` — 9 tautological raises **[R]/[M]**

See §4. `with pytest.raises(X): raise err` where `err` is an `X` constructed two
lines earlier. **[M]** Zero of the 14 guard-disabling mutations reddened
anything in this file.

### 9.3 Two white tests documented as "hand-computed" are blind to the weight formula **[M]**

See §6 caveat 3. `test_white_bc_4_point_quadrature_hand_computed` and
`test_white_bc_axis_z_on_product_quadrature` both PASS under M3 (the `|Ω·n|`
factor dropped). Their inputs are hemisphere-constant, so the measured
functional (a normalised weighted average of a constant) has the entire weight
formula inside its invariance group. **The docstring "explicit hand
calculation" should be corrected** — there is no hand-computed number in the
test; the expected `1.0` follows from normalisation alone.

### 9.4 Three permanently-inert "sentinel" tests behind `except Exception → skip` **[R]/[M]**

**[M]** The boundary harness baseline is `303 passed, **3 skipped**`; **[M]** all
three skips come from one site, with the reason
`2-D mesh construction not available here: tuple index out of range`.

**[R]** `tests/sn/operators/test_bc_extraction_matvec.py:436-445`:

```python
try:
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    if sn_mesh.spatial_shape[1] <= 1:
        pytest.skip("2-D Cartesian vacuum bit-identity needs an ny>1 mesh; …")
except Exception as exc:  # pragma: no cover - construction guard
    pytest.skip(f"2-D mesh construction not available here: {exc}")
```

The `mesh` two lines above is a **1-D** `Mesh1D`, so `spatial_shape` is a
1-tuple and `spatial_shape[1]` raises `IndexError` — the `except Exception`
swallows it. **[R]** The test's own docstring calls it *"the SENTINEL that the
1-D carve did not accidentally perturb the 2-D path"*. It has never run.

Two defects, not one: (a) the sentinel is inert (an acknowledged
"O.4a.2 implementer task" that was committed as a placeholder), and (b) the
bare `except Exception → pytest.skip` means **any future construction bug is
converted into a green skip**, permanently. The intended explicit-skip branch
is unreachable.

### 9.5 `test_bc_equivalence_snapshot.py` is marked `l1` but is a self-generated regression snapshot **[R]/[M]**

**[R]** `pytestmark = [pytest.mark.l1, pytest.mark.regression]` (`:75`) and the
docstring claims *"`@pytest.mark.l1` — cross-implementation L1 (the realizer
path is a structurally-independent re-implementation; the snapshot pins the
numerical reference)"*.

**[R]** But the file's own header records that the cross-implementation half was
**removed**: *"Issue #186 (B3 + β2) removed the `apply` method from every
concrete BC descriptor, eliminating the legacy half. Only the realiser-path
assertion remains; **the snapshots themselves now record realiser-path
outputs**."* A snapshot regenerated from the code it gates is not a
cross-implementation reference at any level — it is exactly what the second
marker says, a `regression` drift gate.

**[M]** The `l1` marker is live: `pytest -m l1 --collect-only` collects all 8.
Under `vv-principles`, this is level conflation (an L1 label on evidence that
produces zero correctness information). The snapshots are genuinely valuable —
**[M]** they reddened under 9 of the 12 leaf mutations in §3.1, the widest net
in the subsystem — but they should carry `regression` alone.

**Adjacent:** `_load_or_skip` **skips rather than fails** when a snapshot file
is missing (`:83-95`). **[M]** All 8 `.npz` are currently git-tracked and all 8
tests pass, so the gate is live today — but a deleted/renamed snapshot silently
disables it rather than reddening. Given §9.4 in the same subsystem, tighten
this to a hard failure.

### 9.6 `tests/geometry/test_boundary_trace_law.py` pins `geometry_map` / `response_kernel` as `None` **[R]**

**[R]** `:95-101` — `test_geometry_map_default_is_none` and
`test_response_kernel_default_is_none` assert `_StubLaw().geometry_map is None`
/ `.response_kernel is None`. These are the ONLY tests touching the two
declared-but-empty properties from review-plan §3, and they pin the *empty*
state. Not a defect (the ABC default IS `None` and testing a default is fine),
but the plan should know: **populating `geometry_map`/`response_kernel` on the
concrete laws is test-neutral** — nothing will red, because nothing asserts
anything about a populated value anywhere. The 14 tests in this file are all
`@pytest.mark.foundation` structural smoke (`hasattr` / `isinstance` / `is None`),
which is correct for an ABC-contract file.

---

## 10. What the 417 bare asserts actually assert [M]

AST classification (`$CLAUDE_JOB_DIR/tmp/classify_asserts.py`), 416 parsed:

| class | count | share |
|---|---:|---:|
| `structural` (`isinstance` / `hasattr` / `is` / `type()`) | 186 | 44.7 % |
| **`numeric`** (the only class that can pin a VALUE) | **121** | **29.1 %** |
| `tag_equality` (`== "reflective"`, the stringly-typed checks) | 46 | 11.1 % |
| `shape_len` (`.shape` / `len()` / `.ndim`) | 36 | 8.7 % |
| `predicate` (a bare bool attr, e.g. `op.is_adjointable`) | 19 | 4.6 % |
| `membership` (`in` / `not in`) | 8 | 1.9 % |

**Value-pinning bare asserts: 121/416 = 29.1 %. Non-value: 295/416 = 70.9 %.**

Adding the function-call assertions (130 `np.testing.*` + 39 `assert_*`
invariant calls, essentially all value/contract gates) the subsystem's overall
value-gating share is **(121 + 130 + 39) / 676 = 42.9 %**.

The 46 `tag_equality` asserts are the test-side shadow of review-plan §3.1's
stringly-typed-dispatch finding — `test_boundary_conditions.py` alone spends 11
of its 17 bare asserts on `sn.bc["xmin"] == "reflective"`-style checks. If the
plan replaces the string dispatch with typed `geometry_map` / `response_kernel`
specifications, those 46 assertions are the migration surface.

---

## 11. Verdict

| brief question | answer |
|---|---|
| 1. Mode-8 `-O` exposure | **REFUTED. 0.0 % of the suite's 676 assertions are inert.** Measured synthetically and on two real boundary files. The only inert bare assert is 1 in a non-collected generator script. |
| 2. Which files exist | 21 enumerated (§2); 390 collected items in the machinery-under-test set, 303 + 3 skipped in the 18-file mutation harness. |
| 3. GATED vs exercised | Every law × method cell is action-gated or refusal-gated, **all mutation-verified** (§3, §5.2). Only SN-periodic is smoke-only, because its action IS the identity. |
| 4. Mutation-verify the strongest claim | 12 leaf mutations + 14 guard-disablings + 5 adjoint mutations = **31 mutations, 30 caught**. The one that passed silently is §9.1 (`catches("ERR-052")`). |
| 5. White BC conservation | **Tested and mutation-live** (§6) — with three honest-scope caveats: procedural not structural independence, operator-level only (white is unreachable from an SN `BC` tag), and two "hand-computed" tests that are blind to the weight formula. |
| 6. Rank-N composition | Gates structure **and** values, both mutation-live (7 / 5 reds). The Marshak value reference is independent of the **walker**, not of the **leaves** — correct as a factorisation, but not "checked against something independent". Rank-N has **no production consumer**. |
| 7. Adjoint coverage | Gated, incl. G-metric reciprocity on 4 geometries. Both refusals gated by exactly one test each. One ungated defensive `MissingAdjoint` raise. |

**Headline for the plan:** the boundary suite is in materially better shape than
the Mode-8 hypothesis suggested — it is the **most mutation-live subsystem I
have measured in this repo**: 30/31 deliberate defects caught, every production
guard live. Its weaknesses are *labelling*, not *coverage*: one stale
`catches` marker, one tautological error-shape file, two over-claiming
docstrings, one level-conflated `l1` marker, and three permanently-skipped
sentinels behind a swallowed exception.

### Recommended plan items (all in-session-fixable, none require new physics)

1. **`catches("ERR-052")`** — re-verify or relocate (§9.1). *Highest priority:
   it is a false coverage edge in the ERR audit.*
2. **`except Exception → pytest.skip`** at
   `tests/sn/operators/test_bc_extraction_matvec.py:445` — either wire the 2-D
   builder the docstring names, or delete the 3 placeholder rows. Never leave a
   bare-`except`-to-skip in a sentinel (§9.4).
3. **`test_bc_equivalence_snapshot.py`** — drop the `l1` marker, keep
   `regression`; make `_load_or_skip` **fail** on a missing snapshot (§9.5).
4. **`test_bc_errors.py`** — delete the 9 tautological
   `with pytest.raises(X): raise err` legs (§4).
5. **Docstring corrections** — `test_white_bc_4_point_quadrature_hand_computed`
   is not a hand calculation (§9.3); `AngularAverageOperator`'s conservation
   claim should note it is verified at the primitive only, never end-to-end
   (§6 caveat 2).
6. **If the plan materialises `geometry_map` / `response_kernel`** — note that
   nothing currently asserts a populated value (§9.6), and that the 46
   `tag_equality` assertions (§10) are the migration surface.
7. **If the plan changes `_reflect_trace`'s outflow discard** — four
   assertions actively DEMAND the outflow rows be zero (§3.1, M8). That is a
   committed contract, not an oversight.
