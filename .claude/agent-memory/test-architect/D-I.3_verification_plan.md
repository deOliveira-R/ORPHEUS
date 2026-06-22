# D-I.3 verification plan — StreamingOperator bare-ndarray retirement

**Wave**: D-I.3 of Depth B (Field-on-FunctionSpace) refactor.
**Surgical target**: `orpheus/sn/operator.py` L1539–L1654 — the bare-`np.ndarray`
arms of `StreamingOperator.apply`. Both the 1-D B1''-packed path and the
2-D `ravel(order='F')`/`solution_to_angular_flux` wrap-and-repack path
retire. After surgery, `StreamingOperator.apply` consumes only
`TimedFullField` and returns only `TimedFullField`.

**Sibling wave state** (read from dispatch brief):
- D-I.1 (c97897d): `CollisionOperator.apply / .solve / .apply_transpose`
  TimedFullField-only. Decomposition gate
  (`tests/sn/test_streaming_operator_decomposition.py`) already migrated.
- D-I.2 (8a8ddbf): `ScatteringOperator.apply(np.ndarray)` retired.
- D-K.1 (400ca33): `SNSolver.L = (StreamingOperator + CollisionOperator)
  = InvertibleOperator`.
- D-K.2 (dadf4e8): `SNStreamingOperator` retired.

**Lesson scaffolding**: this plan invokes L12 (paste-back evidence),
L13 (canonical-type naming), L17 (convention crosswalk), L18 (Pattern 7
producer-side normalisation).

**Augmentation note**: this plan assumes the explorer's dependency audit
returns READY (no production callers still passing bare ndarrays). If
the audit returns BLOCKED, Section 2's deletion verdict for
`test_b1pp_verification.py::test_b1pp_decode_encode_roundtrip` is the
only piece that survives unchanged; the rest of the plan defers behind
the production-caller migration. The risk register (Section 9) treats
this contingency.

---

## Section 1 — Pre-surgery green-gate baseline

Per L12: every claim post-surgery must paste-back against a pre-surgery
stdout snapshot captured BEFORE D-I.3c is dispatched. Three baseline
invocations, run sequentially from the worktree root.

### 1a. Targeted gate (the tests most likely to touch)

```
python -O -m pytest \
  tests/sn/test_streaming_operator.py \
  tests/sn/test_streaming_operator_decomposition.py \
  tests/sn/test_invertible_operator.py \
  tests/sn/test_b1pp_verification.py \
  tests/sn/test_phase_c_gates.py \
  tests/sn/test_2d_l2_matvec_correctness.py \
  tests/sn/test_2d_l2_face_view_unit_source.py \
  -q
```

Capture this stdout as `D-I.3_baseline_targeted.txt`. Every test count
(passed/failed/xfailed/xpassed/skipped) and ULP-sensitive numeric
assertion message becomes the paste-back anchor.

### 1b. Broader SN gate (catches inadvertent breakage)

```
python -O -m pytest tests/sn/ -q \
  --ignore=tests/sn/test_l1_standoff_slab_cylinder.py \
  --deselect tests/sn/test_curvilinear_convergence.py
```

(The L1 slow convergence files are excluded to keep the gate runnable
in a typical session — but see 1c for the MMS convergence pin that
MUST be re-run separately.)

Capture as `D-I.3_baseline_broader.txt`.

### 1c. L1 MMS / convergence pins (load-bearing)

Per L1 pins identified in Section 4 below:
- `tests/sn/test_2d_l2_matvec_correctness.py::test_solve_sn_2d_krylov_homogeneous_reflective_recovers_kinf`
  — pillar: closed-form `k_inf = νΣ_f/Σ_a`. Already in 1a's gate.
- `tests/sn/test_2d_l2_matvec_correctness.py::test_2d_reflective_xy_keff_matches_1d_slab_reflective_analog`
  — pillar: structurally-independent (1-D slab analytic reference).
  Already in 1a's gate.
- Curvilinear MMS convergence (if present):

  ```
  python -O -m pytest tests/sn/test_curvilinear_convergence.py -q
  ```

  Capture stdout as `D-I.3_baseline_l1_mms.txt`. If the file does not
  exist post-Glob, this paste-back is skipped (no pin → no obligation).

---

## Section 2 — Test categorization

Per the dispatch brief deletion criterion:
- **Legacy-contract-only (DELETE)**: tests that exercise the BARE-NDARRAY
  contract via `solution_to_angular_flux*` / `pack_with_traces` for
  reasons unrelated to a mathematical claim production code depends on.
- **Behaviour-pin (MIGRATE)**: tests pinning Resolution A, B1''
  algebraic consistency, GMRES convergence, flat-flux collapse, etc.
  using bare ndarrays as a calling convention. Migrate to TimedFullField.
- **Already-typed (NO-OP)**: tests already TimedFullField; verify they
  still pass.

### `tests/sn/test_streaming_operator.py` (532 LoC)

| Class / test | Verdict | Note |
| --- | --- | --- |
| `TestCapabilities.*` (L161–L200) | MIGRATE | Capability advertising. Doesn't call `apply`; needs no ndarray. Already type-clean — no change. |
| `TestConstruction.test_construct_with_sn_mesh_and_sigma_t` (L212) | NO-OP | No apply call. |
| `TestApplyShape.test_apply_preserves_shape` (L229) | **DELETE** | Pure ndarray-shape contract via `_packed_psi` (L233). The TimedFullField shape contract is already pinned by `TestCompositeInvariants.test_returns_timed_full_field` (L385). |
| `TestApplyShape.test_apply_returns_packed_vector` (L238) | **DELETE** | Same — asserts `out.ndim == 1` on the retired packed path. Replaced structurally by `TestCompositeInvariants.test_returns_timed_full_field` (TimedFullField is the new shape contract). |
| `TestLinearity.test_apply_zero_returns_zero` (L274) | MIGRATE | Linearity pin. Math claim survives; rewrite using `_random_composite` (L359) or `sn_mesh.zeros_timed_full_field()`. |
| `TestLinearity.test_apply_is_linear` (L284) | MIGRATE | Linearity pin. Same migration: `alpha*state1 + beta*state2` via TimedFullField arithmetic (already supported per `test_LC_apply_composite` L510). |
| `TestSumCapabilities.*` (L317–L338) | NO-OP | No apply call. |
| `TestCompositeInvariants.*` (L384–L496) | NO-OP | Already TimedFullField. Pins post-surgery contract verbatim. |
| `TestOperatorAlgebraCompositionUnderTimedFullField.*` (L510+) | NO-OP | Already TimedFullField. |

**Helper `_packed_psi` (L128–L145)**: DELETE after the two TestLinearity
tests migrate. The helper exists solely to size against the retired
packed contract.

### `tests/sn/test_streaming_operator_decomposition.py` (333 LoC)

Per D-I.1 commit message (the file already migrated to TimedFullField),
`grep -n "TimedFullField\|np.ndarray"` confirms: all three tests
(`test_bit_exact_uniform_sigma_t` L147, `test_L_apply_equals_subtractive_form`
L227, `test_subtractive_L_differs_from_matvec_at_zero_sigma_t` L294)
construct `TimedFullField` explicitly at fixture time (L157, L234, L303).

Verdict: **ALL NO-OP**. The decomposition gate's Resolution A bit-exact
contract (Section 3) re-verifies verbatim post-surgery.

### `tests/sn/test_invertible_operator.py` (926 LoC)

`grep` confirms: TimedFullField imported at L55, helpers `_random_state`
(L61) and `_constant_state` (L82) return TimedFullField. The legacy
contract is present only as a NEGATIVE assertion:
`test_solve_rejects_bare_ndarray` (L310) — pins that the typed solver
REJECTS bare ndarrays. This contract STAYS post-surgery (the apply
narrows similarly; the rejection contract holds equally).

Verdict: **ALL NO-OP** except `test_solve_rejects_bare_ndarray` which
STRENGTHENS (it now also applies to `.apply`, not just `.solve`). Add
a sibling test:

```python
def test_apply_rejects_bare_ndarray(self) -> None:
    """D-I.3 — StreamingOperator.apply must raise on bare ndarrays."""
    sn = _slab(nx=4)
    sigma_t = _sig_t(sn, ng=2)
    L = StreamingOperator(sn, sigma_t)
    with pytest.raises(TypeError, match="TimedFullField"):
        L.apply(np.zeros(L.n_unknowns))
```

This is the typed-contract NEGATIVE test (vv-principles anti-pattern
#11: contract-validation methods need positive + negative tests).

### `tests/sn/test_b1pp_verification.py` (302 LoC) — the densest legacy site

| Test | Verdict | Note |
| --- | --- | --- |
| `test_b1pp_lplusc_is_full_rank` (L116) | **MIGRATE** | Pins MATH: (L+C) full-rank under B1''. Currently builds a dense matrix by probing with unit basis vectors via `matvec(e)` where `e: np.ndarray`. **Migration**: the dense-matrix-probe is unavoidable for a rank test (need a square matrix); use the **face-flat-buffer view** of `BoundaryFlux` + `AngularFlux.values.ravel()` to construct unit basis vectors *as TimedFullField components*, then `np.linalg.matrix_rank(M)` where `M[:, k]` is the concatenation of bulk-ravel + boundary-ravel of `L.apply(e_k_typed) + C.apply(e_k_typed)`. The dimensionality changes (typed B1'' may have a different layout from bare-packed B1''); the math claim (full-rank) is layout-agnostic. |
| `test_b1pp_constant_flux_collapses_to_collision` (L167) | **MIGRATE** | Pins MATH: flat ψ ⇒ (L+C)·1 = σ_t·1 at cell slots, 0 at face slots. **Migration**: build `state = sn_mesh.ones_timed_full_field()` (or construct explicitly with `state.bulk.values[...] = 1.0` and `state.boundary.values[...] = 1.0`), then assert `np.allclose(out.bulk.values, sigma_t_val)` and `np.allclose(out.boundary.values, 0.0)` (where the boundary IS the face-residual buffer per D-G). The packed `out[:n_cell_scalars]` / `out[n_cell_scalars:]` split disappears — replaced by `out.bulk` / `out.boundary` (Pattern 4: illegal-states-unrepresentable — the slot distinction becomes a type distinction, not an array-slicing convention). |
| `test_b1pp_lplusc_gmres_converges_fp_noise` (L224) | **MIGRATE** | L1 GMRES convergence pin. **Critical**: GMRES requires a `scipy.sparse.linalg.LinearOperator` over **a flat ndarray**. The migration MUST preserve this — wrap TimedFullField ↔ flat-ndarray at the GMRES boundary via `state.flat_view()` / `TimedFullField.from_flat(flat, mesh, history_depth=...)` (assuming D-G shipped these; if not, this test is BLOCKED until D-G ships the round-trip helper that the rest of D-I assumes). Alternative: pin GMRES on **`InvertibleOperator.solve`** (`A.solve(rhs_state)`) since `A = L + C` already provides a typed Krylov dispatch. Prefer this — it's the production path. **Verdict**: rewrite to drive GMRES via `(L+C).solve(q_typed)` and assert convergence info + residual. |
| `test_b1pp_decode_encode_roundtrip` (L279) | **DELETE** | Pins the BARE-NDARRAY contract VIA `solution_to_angular_flux_with_traces` + `pack_with_traces`. These helpers are the machinery being retired (per dispatch brief). The test exists solely to gate the contract that is going away. Per Section 2 of the brief: tests pinning the bare-ndarray contract via the now-retired `solution_to_angular_flux*` / `pack_with_traces` family are by definition pinning machinery that is going away. DELETE alongside the surgery. |

### `tests/sn/test_phase_c_gates.py` (889 LoC)

Per the D-K.5 migration comment at L21–L25 ("Migrated from the retiring
:class:`SNStreamingOperator` (packed-vector matvec API) to the composite
algebra :class:`StreamingOperator` + ... :class:`TimedFullField`"),
this file is **ALL NO-OP**. All tests construct via TimedFullField
helpers (`_build_timed_full_field` L146, `_flat_psi_timed_full_field` L242).

The xfail at L353 (`test_apply_apply_transpose_reciprocity_under_sweep_frame`)
is a pre-existing reciprocity contract waiting on Wave J; D-I.3 does
NOT touch it.

### `tests/sn/test_2d_l2_matvec_correctness.py` (305 LoC) and `test_2d_l2_face_view_unit_source.py` (239 LoC)

Both files already use `mesh.zeros_timed_full_field()` (Section 2.2:
`test_apply_vs_sweep_2d_residual_cancellation` L171–L207 explicitly
builds a TimedFullField at L183). **ALL NO-OP**.

The deferred Test 3.1 MMS L1 pin (Issue #210 — `_apply_2d_cartesian`
vectorization) remains deferred; D-I.3 does NOT unblock it.

### Summary table — deletion verdict

| File | DELETE | MIGRATE | NO-OP | NEW |
| --- | --- | --- | --- | --- |
| `test_streaming_operator.py` | 2 (TestApplyShape) + `_packed_psi` helper | 2 (TestLinearity) | rest | 0 |
| `test_streaming_operator_decomposition.py` | 0 | 0 | 3 | 0 |
| `test_invertible_operator.py` | 0 | 0 | all | 1 (apply rejects bare ndarray) |
| `test_b1pp_verification.py` | 1 (roundtrip) | 3 | 0 | 0 |
| `test_phase_c_gates.py` | 0 | 0 | all | 0 |
| `test_2d_l2_matvec_correctness.py` | 0 | 0 | all | 0 |
| `test_2d_l2_face_view_unit_source.py` | 0 | 0 | all | 0 |

Net: **3 deletions** (2 TestApplyShape + 1 b1pp roundtrip) + **5
migrations** (2 TestLinearity + 3 b1pp behaviour-pins) + **1 new
negative test** (apply rejects bare ndarray).

---

## Section 3 — Resolution A bit-exact decomposition gate

The load-bearing correctness gate for the operator algebra is
`tests/sn/test_streaming_operator_decomposition.py::TestResolutionABitExact`
(per D-I.1, gates `(L+C).apply ≡ M` bit-exactly on bulk AND boundary).

Per Section 2 above, the file is ALREADY TimedFullField (D-I.1
shipped the migration). After D-I.3:
- No source-of-the-test changes (the test fixture builds `TimedFullField`
  inline; the apply path narrows from "consume TimedFullField OR
  ndarray" to "consume TimedFullField only" — no observable change).
- The gate re-runs verbatim.

Post-surgery verification command:

```
python -O -m pytest tests/sn/test_streaming_operator_decomposition.py -v
```

Expected: 3 tests pass, no diff from baseline. Paste-back this against
1a's baseline section for the decomposition file.

Per `vv-principles` §Bit-identity vs principled-equivalence: D-I.3
does NOT change the FP reduction tree (the typed path was already the
post-decode path inside the retired branch; the surgery only removes
the decode-then-call indirection, not the math). Expect bit-identity.
If bit-identity breaks, investigate — the surgery has changed
something other than what the dispatch brief states.

---

## Section 4 — L1 MMS convergence pins

Identified L1 MMS / closed-form pins on `StreamingOperator.apply`:

| Test | File:line | Pillar | Post-surgery status |
| --- | --- | --- | --- |
| `test_solve_sn_2d_krylov_homogeneous_reflective_recovers_kinf` | `tests/sn/test_2d_l2_matvec_correctness.py:137` | Closed-form `k_inf = νΣ_f/Σ_a = 1.875` | Already typed. Re-runs verbatim. |
| `test_2d_reflective_xy_keff_matches_1d_slab_reflective_analog` | `tests/sn/test_2d_l2_matvec_correctness.py:214` | Structurally-independent (1-D slab analytic ref) | Already typed. Re-runs verbatim. |
| `test_apply_vs_sweep_2d_residual_cancellation` | `tests/sn/test_2d_l2_matvec_correctness.py:171` | Algebraic identity (`OperatorSum.apply` distributes) | Already typed; `@catches("ERR-026")`. Re-runs verbatim. |
| `test_b1pp_lplusc_gmres_converges_fp_noise` | `tests/sn/test_b1pp_verification.py:224` | B1'' algebraic consistency (GMRES floor) | MIGRATE (Section 2) — preserves the convergence claim through `(L+C).solve` |

**Issue #210 status (noted, not blocking)**: Test 3.1 MMS L1
convergence pin on `_apply_2d_cartesian_l2` is deferred pending
vectorization. D-I.3 does NOT block on this; the pin can land
independently after `_apply_2d_cartesian` is vectorized.

**Curvilinear MMS convergence** (slab/sphere/cylinder): if
`tests/sn/test_curvilinear_convergence.py` exists, run it as part of
1c's baseline and re-run post-surgery for paste-back. If the file is
absent, no pin is in scope here (the curvilinear convergence claim is
verified elsewhere, e.g. by L1 standoff cylinder/slab vs analytic).

Post-surgery L1 gate command (combines all L1 pins on touched paths):

```
python -O -m pytest \
  tests/sn/test_2d_l2_matvec_correctness.py \
  tests/sn/test_b1pp_verification.py \
  tests/sn/test_l1_standoff_slab_cylinder.py \
  -v -m "l1 or foundation"
```

Paste-back against `D-I.3_baseline_l1_mms.txt`.

---

## Section 5 — Production-path L0 anchor

The L0 streaming-equilibrium analytical anchor (`φ = Q/Σ_t(1-c)`,
`ψ_n = φ/Σw`) is the universal probe for streaming + collision +
scattering (Signature 4 in `numerical-bug-signatures`). Per dispatch
brief: the SNSolver pipeline routes through `InvertibleOperator.solve`
/ `.apply` consuming TimedFullField — no production-path change at
D-I.3, only API narrowing.

Tests pinning this anchor (search target: `Q/sigma_t`,
`streaming.*equilibrium`, fixed-source flat-flux):
- `tests/sn/test_streaming_equilibrium.py` (likely path — verify via
  `Glob` at execution time).
- `tests/sn/test_quadrature.py::TestL0TermVerification::test_per_ordinate_flat_flux_consistency`
  (referenced in `numerical-bug-signatures` Signature 1).

These tests construct via the production `solve_sn` API, which already
consumes TimedFullField post-D-K.1. They re-run verbatim. Verify
by inclusion in the broader gate (1b).

---

## Section 6 — Cross-geometry coverage matrix

Coverage matrix post-surgery — confirms the carve preserves every
geometry through a typed-path test plus a multi-group, heterogeneous
case (per `vv-principles` anti-patterns #3, #4 and L2: 1-group proves
nothing).

| Geometry | Carrier | File:test | ≥2G? | Heterogeneous? |
| --- | --- | --- | --- | --- |
| 1-D slab (cartesian, ny=1, curv=None) | TimedFullField | `tests/sn/test_streaming_operator.py:385` (`test_returns_timed_full_field[slab]`) | 2G via `_sig_t_uniform(ng=2)` | NO — but covered by `test_streaming_operator_decomposition.py:227` (`test_L_apply_equals_subtractive_form[slab]`) using `_random_state` + heterogeneous σ_t |
| 1-D slab (algebraic identity) | TimedFullField (packed→typed via decomposition gate) | `test_streaming_operator_decomposition.py:147` (`test_bit_exact_uniform_sigma_t[slab]`) | 2G | Uniform σ_t (foundation pin); heterogeneous covered by L1 `test_l1_standoff_slab_cylinder.py` |
| 1-D sphere (curvilinear) | TimedFullField | `tests/sn/test_streaming_operator.py:385` (`test_returns_timed_full_field[sphere]`) | 2G | NO at L0 |
| 1-D sphere (algebraic) | TimedFullField | `test_streaming_operator_decomposition.py:147` (`test_bit_exact_uniform_sigma_t[sphere]`) | 2G | Uniform; heterogeneous via L1 cylinder/sphere standoff |
| 1-D sphere (curvilinear correctness) | TimedFullField | `test_phase_c_gates.py:258` (`test_apply_curvilinear_per_ordinate_flat_flux_residual`) — `@verifies("dd-curvilinear-scalar") @catches("ERR-026")` | 2G via `placeholder_materials(ng=2)` | YES (parametrised σ_t ∈ {0, 0.5}) |
| 1-D cylinder | TimedFullField | `test_streaming_operator.py:385` (`test_returns_timed_full_field[cylinder]`) + `test_streaming_operator_decomposition.py:147` (`[cylinder]`) | 2G | Cylinder heterogeneous: covered by L1 cylinder standoff |
| 2-D Cartesian | TimedFullField | `test_2d_l2_matvec_correctness.py:171` (`test_apply_vs_sweep_2d_residual_cancellation`) — random TimedFullField state + non-zero scatter + 2g mixture | 2G via `get_mixture("A", "2g")` | Vacuum with non-uniform random state activates redistribution |
| 2-D Cartesian (k_inf) | TimedFullField (via `solve_sn`) | `test_2d_l2_matvec_correctness.py:137` | 2G | Homogeneous (closed-form anchor) — heterogeneous covered by `test_2d_reflective_xy_keff_matches_1d_slab_reflective_analog` L214 |

**Coverage gaps surfaced**:
- 1-D slab heterogeneous + L0 typed: not directly tested at TimedFullField level. ACCEPTABLE — the heterogeneous claim is reachable through L1 standoff against analytic. If post-surgery L1 standoff drifts, this is the gap to widen.
- 1-D sphere/cylinder heterogeneous + ≥2G + typed at L0: covered transitively by `test_phase_c_gates.py:258`. ACCEPTABLE.

No new coverage tests required by D-I.3 (the surgery does not introduce
new mathematical territory; it only narrows the input contract).

---

## Section 7 — Convention crosswalk (per L17)

For the typed↔packed bridge being collapsed at D-I.3:

| Subsystem | Input convention | Internal convention | Output convention |
| --- | --- | --- | --- |
| `StreamingOperator.apply` **pre-D-I.3 1-D bare-ndarray** | packed flat: `[ψ_cell.T.ravel(F), ψ_face_outer.T.ravel(F), ψ_face_inner.T.ravel(F)]` | wrap → `TimedFullField` internally at L1603–L1622 via `solution_to_angular_flux_with_traces` decoder; call `_apply_timed_full_field`; gather face residuals into m_face_outer/m_face_inner buffers | pack → flat via `pack_with_traces` (L1644), then subtract `σ_t · ψ_packed[:n_cell_scalars]` at cell-only slots (L1647–L1654) |
| `StreamingOperator.apply` **pre-D-I.3 2-D bare-ndarray** | packed flat: `ψ.ravel(order='F')` per `(ng, nx*ny, n_ord)` | wrap via `solution_to_angular_flux` (L1571) + zero `BoundaryFlux`; call `_apply_timed_full_field` | repack via `out_typed.bulk.values.ravel(order='F')` (L1586–L1588) — note: 2-D path returns plain ndarray (no σ_t·ψ subtraction at packed level; Resolution A subtraction already happens inside `_apply_2d_cartesian`) |
| `StreamingOperator.apply` **post-D-I.3** | `TimedFullField` only | `TimedFullField` throughout | `TimedFullField` only |
| `InvertibleOperator.apply` (consumer) | `TimedFullField` (per D-I.1) | delegates: `L.apply(state) + C.apply(state)` (per `OperatorSum.apply`) | `TimedFullField` |
| `SourceIteration.solve` / `solve_sn` (production) | `TimedFullField` (per D-H.1c, D-K.1) | n/a (calls through `(L+C).solve`) | `TimedFullField` |
| Bare-ndarray test callers (pre-surgery) | packed flat | wraps to TimedFullField at the typed adapter inside `apply` | packed flat |
| Bare-ndarray test callers (post-surgery, post-migration per Section 2) | **REMOVED** — every test consumes `TimedFullField` | — | — |

**Crosswalk verdict** (per L18 Pattern 7):
- The convention-bridge collapse REMOVES the pack/unpack code paths.
  The producer (`_apply_timed_full_field`) is the SINGLE convention
  fix-point; consumers narrow to a single input type.
- No new producer-consumer mismatch is created. The pattern's
  invocation is correctly oriented: convention lives at the producer,
  consumers see only the post-normalised type.
- **Specifically**: the σ_t · ψ subtraction was DUPLICATED between the
  packed path (L1647 — subtracts at packed cell slots only) and the
  typed path (`_apply_timed_full_field` — subtracts at cell-centres
  per Resolution A). Removing the packed path eliminates the
  duplication; only the typed-path subtraction remains. Per Pattern 2
  (single source of truth), this is a beneficial collapse. The
  bit-exact decomposition gate (Section 3) confirms the typed-path
  subtraction is the principled formulation.

**Risk surface from the crosswalk**: NONE. Every convention narrows;
no new convention is introduced.

---

## Section 8 — Verification protocol order

The exact ordering for D-I.3 verification, with paste-back gates:

1. **Pre-surgery baseline** (Section 1a + 1b + 1c).
   - Owner: main agent (NOT a sub-agent — paste-back must be captured
     verbatim with no summarisation per L12).
   - Output: three stdout snapshots saved under
     `.claude/agent-memory/test-architect/D-I.3_baseline_*.txt`.

2. **D-I.3c sub-agent dispatch — test migration**.
   - **Owner**: dispatched method-implementer OR direct main-agent
     (Section 9 risk b favours main-agent direct; see
     `feedback_no_method_implementer_for_surgical_carves`).
   - **Scope**: apply Section 2's deletions + migrations exactly.
     Three deletions, five migrations, one new negative test.
   - **Paste-back required**: post-migration, re-run command 1a (the
     targeted gate). Tests in Section 2's MIGRATE rows must pass
     against the (still-existing) bare-ndarray apply path. Tests in
     Section 2's DELETE rows must be GONE from the test count. The
     new negative test must PASS (because the bare-ndarray rejection
     is the post-surgery behaviour — wait, that's wrong; pre-D-I.3d
     the apply still accepts bare ndarrays. The new test must
     therefore be marked `xfail(strict=True)` until D-I.3d ships).
   - Per the `feedback_vv_tagging.md` convention: the new
     `test_apply_rejects_bare_ndarray` ships at D-I.3c with
     `@pytest.mark.xfail(reason="D-I.3d narrows the API", strict=False)`,
     flips to passing at D-I.3d when the bare-ndarray arm is removed.

3. **Main-agent independent re-run of D-I.3c paste-back** (per L12).
   - Owner: main agent.
   - Action: re-execute command 1a; diff against the sub-agent's
     paste-back. Any mismatch ⇒ investigate before D-I.3d.

4. **D-I.3d operator surgery** — remove L1539–L1654 of
   `orpheus/sn/operator.py`.
   - Owner: main-agent DIRECT (per
     `feedback_no_method_implementer_for_surgical_carves`: surgical
     carves on operator algebra are too important to lose real-time
     correction).
   - Action: delete the `else` branch at L1539. Confirm the
     dispatch comment + xfail flag at the new negative test flip to
     `strict=True`.

5. **Post-surgery full gate**:
   - Re-run 1a, 1b, 1c. Paste-back against baselines from step 1.
   - Re-run Resolution A decomposition gate
     (`tests/sn/test_streaming_operator_decomposition.py` per Section 3).
   - Re-run L1 MMS gate per Section 4.
   - Verify the new negative test now passes (strict=True flip).
   - Owner: main agent. Paste-back filed as
     `D-I.3_postsurgery_*.txt`.

6. **D-J cascade orphan check** (test-architect gate, CAN BE DEFERRED).
   - Action: dispatch test-architect (this agent) to audit whether
     `solution_to_angular_flux`, `solution_to_angular_flux_with_traces`,
     `pack_with_traces`, `EquationMap` retain any callers outside
     `orpheus/sn/operator.py`. If yes — those callers need their own
     retirement waves (D-J.1, D-J.2, …). If no — the helpers
     themselves retire at D-J.0.
   - Owner: dispatched test-architect or general-purpose with explorer
     dependency audit.

---

## Section 9 — Risk surface

Top risks, ranked by impact × likelihood:

### (a) Hidden production caller passes bare ndarrays — BLOCKING

**Mechanism**: a production caller in `orpheus/sn/iteration.py`,
`solver.py`, or a Krylov adapter still passes a bare ndarray through
`StreamingOperator.apply` (e.g. via a `scipy.sparse.linalg.LinearOperator`
matvec wrapper). The surgery breaks production.

**Mitigation**: Explorer dependency audit (running in parallel per
the dispatch brief). This plan assumes the audit returns READY.
If the audit returns BLOCKED, D-I.3d does not ship until the offending
production caller migrates to `TimedFullField` (or to `(L+C).solve`
which already accepts TimedFullField).

**Test signal**: the broader gate (1b) is the catch-net. If a
production path silently coerces to ndarray, the broader gate fails
post-surgery despite the targeted gate passing.

### (b) Test-deletion criterion ambiguity — MEDIUM, mitigated

**Mechanism**: a test pinning math via bare-ndarray calling convention
gets deleted instead of migrated. Math coverage shrinks invisibly.

**Mitigation**: Section 2 explicitly categorises every test in the
seven-file gate. The verdict table makes the deletion criterion
adjudicable per-test. No test gets deleted without a named MIGRATE
counterpart, except the b1pp roundtrip which is unambiguous (it
pins the contract being retired, not a math claim).

**Audit signal**: the post-surgery test count must equal pre-surgery
count minus 3 (deletions) plus 1 (new negative test) = baseline − 2.
Any other delta signals an undocumented deletion.

### (c) FP-non-associativity drift on Resolution A bit-exact gate — LOW

**Mechanism**: the bare-ndarray path performs σ_t · ψ subtraction at
the packed cell-slice level (a flat indexed array). The typed path
performs σ_t ⊙ ψ.bulk subtraction as a `(ng, nx, ny, N)` broadcast
(per `_apply_timed_full_field`). Floating-point reduction tree
differs; bit-identity may fail at ULP level.

**Mitigation**: this drift would manifest in `test_streaming_operator_decomposition.py::test_bit_exact_uniform_sigma_t` (Section 3 gate)
which already runs the TYPED path bit-exactly. Pre-D-I.3 that test
passed (per D-I.1). Post-D-I.3 the test exercises the same typed path
— no FP tree change. Bit-identity is preserved.

**Backstop**: if bit-identity does break unexpectedly, the three
criteria of `vv-principles` §Bit-identity vs principled-equivalence
apply: (1) the new formulation IS principled (per-cell σ_t·ψ
subtraction is the named Resolution A operation); (2) the new value
agrees with `k_inf` analytical limit (Section 4 anchor); (3) drift
bounded by `(iteration_count) × ULP`. Relax the regression contract
for the affected test only, with the documented relaxation.

---

## Closing — the canonical types this plan names

Per L13 (sub-agent briefs MUST name canonical types):

- **`TimedFullField`** — the composite carrier consumed AND returned
  by every operator after D-I.3.
- **`AngularFlux.from_mesh(values, mesh)`** — bulk-only factory used
  for fixture construction in migrated tests.
- **`BoundaryFlux.zeros_for_sn_mesh(mesh)`** — zero-boundary factory
  for `_random_composite`-style fixtures (already in `_random_composite`
  at `test_streaming_operator.py:359`).
- **`sn_mesh.zeros_timed_full_field(history_depth=2)`** — the
  one-liner for "give me a typed zero state with default history"
  used throughout `test_phase_c_gates.py` and the migrated
  `TestLinearity` tests.
- **`_history=()` and `history_depth=2`** — composite-construction
  kwargs preserved by the surgery (Section 3 bit-exact gate verifies
  history pass-through).

These are the types every D-I.3c migration touches; the brief to any
dispatched implementer must name them verbatim.

