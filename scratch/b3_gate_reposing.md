# B3 — re-posing the outflow-discard gates (C-1)

Verification design memo for campaign phase **B3** of the boundary-machinery
review (`.claude/plans/boundary_machinery_review.md` §B3, constraint **C-1**).
Branch `refactor/operator-strategy-layers` @ `acb46245`, repo root
`/Users/rodrigo/git/nuclear/ORPHEUS` (main checkout, not a worktree).
Canonical invocation `.venv/bin/python -O -m pytest`, SERIAL.

Evidence key: **[M]** measured in this session (pytest stdout pasted) ·
**[R]** read (quoted + `file:line`) · **[G]** grep/exhaustive.

**Scope discipline.** Nothing under `orpheus/` or `tests/` was written. Every
mutation is an in-process monkeypatch installed by a pytest plugin under
`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/`; no `git checkout` / `restore` /
`stash` / `clean` was run on any path.

**Instruments built (all under `$CLAUDE_JOB_DIR/tmp/`):**

| file | what it is |
|---|---|
| `b3plugin.py` | the B3 **simulation** — installs the narrowing in-process, plus 7 mutation variants |
| `b3/probe_narrow3.py` | the narrowing-equivalence probe on real `SNMesh` fixtures |
| `b3/probe_strict.py` | does a reduced `PermutationOperator` refuse a full-face input? |
| `b3/probe_variants.py` | sha256 **positive controls** — proves each mutation actually bit |
| `b3/test_reposed.py` | the **proposed** re-posed gates, runnable, so their colour is measured before you land them |
| `b3/run_matrix.sh` | the gate-colour matrix runner |

---

## The change under design, stated precisely

**[R]** `orpheus/sn/operators/boundary.py:351-359` — today, per face:

```python
face_in = boundary.face_view(face)                 # the WHOLE face slot
sel = rows[face] if rows is not None else trace.inflow_indices_for_face(face)
full = law.apply(face_in)                          # law : full face -> full face
out_boundary.face_view(face)[sel] = full[sel]      # <- the outflow rows are discarded
```

**[R]** `orpheus/diffusion/operators.py:611-614` — the honest domain:

```python
for face, law in self.face_laws.items():
    out_boundary.face_view(face)[ScalarTraceSpace.INFLOW_ROW] = (
        law.apply(trace.outflow_view(face))         # law : Gamma_+ -> Gamma_-
    )
```

B3 narrows SN to match: the realized law's **domain becomes `Γ₊`**, its
**codomain `Γ₋`**.

### A structural fact the plan does not state, and every gate below depends on

**[R]** `orpheus/numerics/spaces/angular_trace_space.py:432-448` — the face slot
is a **three-way** split, not two:

```python
def inflow_indices_for_face(self, face):   return np.flatnonzero(row < -TANGENTIAL_EPS)
def outflow_indices_for_face(self, face):  return np.flatnonzero(row > +TANGENTIAL_EPS)
```

Tangential ordinates (`|Ω·n| ≤ ε`) are in **neither** set. **[M]** On the
`cyl_reflective` gate fixture (`Quadrature.product(n_mu=2, n_phi=4)`),
**4 of 8 ordinates are tangential**:

```
### cyl_reflective
  face=xmax  kind=reflective   shape=(8, 1) in=2 out=2 tangent=4 NARROW-BIT-ID=True max|d|=0.000e+00
           tangent rows of the law image: max|.|=1.846e+00  (discarded by the sel write)
```

So `outflow ≠ complement(inflow)`, and **the current gate's `got[outflow] == 0`
leaves the 4 tangent rows completely unguarded** — the law's image there is
non-zero (1.846) today and is silently dropped. §D2/D3 measure the consequence.

---

## D1 — the gates, located and characterised

### D1.1 The plan's count is right about the ASSERTION and wrong about the BLAST RADIUS

**[M]** Re-running the auditor's `reflect_nosel` mutation (M8: drop the `sel`
row restriction, i.e. make `B` preserve the outflow rows) over
`tests/{geometry,sn/operators,diffusion,numerics}`:

```
48 failed, 2174 passed, 5 skipped, 2 xfailed, 34 warnings in 217.21s (0:03:37)
```

The review says "4 tests plus the 2-D balance gate". Both are confirmed, and
both are **one test function each**:

| # | gate | file:line | the assertion |
|---|---|---|---|
| **G1** | `TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row` — **4 params** (`slab_vacuum_reflective`, `slab_reflective_reflective`, `sphere_reflective`, `cyl_reflective`) | `tests/sn/operators/test_sn_boundary_operator.py:123-158` | two legs: `assert_array_equal(got[inflow], bc.apply(face_slot)[inflow])` (:152) **and** `assert not got[outflow].any()` (:155), msg *"B emitted non-zero on the outflow row"* |
| **G2** | `TestBoundaryResidual2DDrivesToZero::test_reflective_boundary_balance_at_convergence` | `tests/sn/operators/test_bc_extraction_2d.py:415-467` | `worst_defect < 10 × solver_tol` where `defect = |ψ.inflow − (B·ψ)[inflow]|`, preceded by the `keff == 1.875` precondition |

**So "4 tests" = 4 parametrisations of ONE function. The real count is 2 test
functions / 5 test items.** That matches C-1's intent; state it as
*5 items, 2 functions* so the re-write is not mis-scoped.

**The correction that matters:** the outflow discard is **asserted** by G1
and **consumed** by 43 more items. The other 43 reds are not contract
assertions — they are value consumers (`solve_sn` reads
`reflect_into_inflow`, the adjoint reciprocity suite reads `B.H`). Under B3
the composite value does **not** change (§D4), so those 43 must stay green;
they are the bit-identity gate's audience, not C-1's rewrite list.

### D1.2 What each of the 5 items catches TODAY (measured)

| mutation (in-process, bite-checked) | G1 | G2 |
|---|---|---|
| M8 `reflect_nosel` — write the law's full-face image | **RED** ×4, *"B emitted non-zero on the outflow row"* | **RED** — `precondition: 2-D reflective keff = … ≠ k_inf 1.875` |
| M1 reflective perm → identity | RED (review's audit) | RED |
| M2 reflective perm → `np.roll(perm,1)` | RED | **GREEN** — a rolled permutation preserves the balance functional (`vv` Mode 12; recorded in the review §7C.2) |

Verbatim G1 red under M8:

```
FAILED tests/sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[slab_vacuum_reflective] - AssertionError: slab_vacuum_reflective face 'xmin': B emitted non-zero on t...
FAILED tests/sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[slab_reflective_reflective] - AssertionError: slab_reflective_reflective face 'xmin': B emitted non-zero ...
FAILED tests/sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[sphere_reflective] - AssertionError: sphere_reflective face 'xmax': B emitted non-zero on the ou...
FAILED tests/sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[cyl_reflective] - AssertionError: cyl_reflective face 'xmax': B emitted non-zero on the outfl...
FAILED tests/sn/operators/test_bc_extraction_2d.py::TestBoundaryResidual2DDrivesToZero::test_reflective_boundary_balance_at_convergence
```

### D1.3 THREE more items C-1 must adopt — the plan's list is incomplete

Under the **honest** full-fidelity B3 (below), three further items break for
the SAME reason G1 does — they feed the realized law a **full-face** array:

| # | item | file:line | why it breaks |
|---|---|---|---|
| **G3** | `TestApplyTransposeCapability::test_adjointable_when_all_faces_support` — **4 params** | `tests/sn/operators/test_sn_boundary_operator.py:188-219` | reference is `bc.apply_transpose(masked)` with `masked = np.zeros_like(face_in)` — a **full-face** array (`:211-214`) |
| **G4** | `TestApplyTransposeCapability::test_adjointability_drops_when_a_face_lacks_it` | same file, `:221-251` | the `_NoTransposeLaw` **duck-typed surrogate** (`:227`) defines `apply(x) = x`, a full-face identity. Invisible to grep and to the graph — this is exactly the **FOURTH search** B2 added to the retirement audit (plan §B2.2, *duck-typed test surrogates*) |
| **G5** | `Test188WiringContracts::{test_curvilinear_realize_boundary_law_routes_through_realizer, test_1d_reflective_faces_realized_as_permutation_tp}` + 6 in `test_snmesh_realizer_wiring.py` | `tests/geometry/test_bound_compat.py`, `tests/sn/operators/test_snmesh_realizer_wiring.py:309…` | they assert the realized object **IS** a `TensorProductOperator` / `PermutationOperator` on the full ordinate set, and probe it with `np.ones(quad.N)` |

**[M] THE DEFINITIVE BLAST RADIUS** — the honest narrowing
(`ORPHEUS_B3=realizer`, strict domain guard) over
`tests/{geometry,sn,diffusion,numerics,transport} -m "not slow"`:

```
16 failed, 4258 passed, 6 skipped, 111 deselected, 57 xfailed, 21 warnings in 989.78s (0:16:29)
```

```
FAILED geometry/test_bound_compat.py::Test188WiringContracts::test_curvilinear_realize_boundary_law_routes_through_realizer
FAILED geometry/test_bound_compat.py::Test188WiringContracts::test_1d_reflective_faces_realized_as_permutation_tp
FAILED sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[slab_vacuum_reflective]
FAILED sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[slab_reflective_reflective]
FAILED sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[sphere_reflective]
FAILED sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[cyl_reflective]
FAILED sn/operators/test_sn_boundary_operator.py::TestApplyTransposeCapability::test_adjointable_when_all_faces_support[slab_vacuum_reflective]
FAILED sn/operators/test_sn_boundary_operator.py::TestApplyTransposeCapability::test_adjointable_when_all_faces_support[slab_reflective_reflective]
FAILED sn/operators/test_sn_boundary_operator.py::TestApplyTransposeCapability::test_adjointable_when_all_faces_support[sphere_reflective]
FAILED sn/operators/test_sn_boundary_operator.py::TestApplyTransposeCapability::test_adjointable_when_all_faces_support[cyl_reflective]
FAILED sn/operators/test_snmesh_realizer_wiring.py::test_2d_cartesian_vacuum_xmin_masks_only_inflow
FAILED sn/operators/test_snmesh_realizer_wiring.py::test_2d_cartesian_reflective_ymax_returns_permutation
FAILED sn/operators/test_snmesh_realizer_wiring.py::test_2d_reflective_y_face_builds_y_axis_permutation
FAILED sn/operators/test_snmesh_realizer_wiring.py::test_1d_cartesian_vacuum_right_masks_only_inflow
FAILED sn/operators/test_snmesh_realizer_wiring.py::test_1d_spherical_vacuum_routes_through_realizer
FAILED sn/operators/test_snmesh_realizer_wiring.py::test_1d_cylindrical_one_boundary_outer_reflective
```

**All 16 are SIGNATURE breaks; ZERO are value breaks.** Every solver,
adjoint-reciprocity, DSA and Gauss-Seidel consumer in `tests/sn/**` stays
green — the strongest confirmation of §D4's bit-identity claim. Earlier
drafts of this memo reported 18; two of those were artefacts of the first
simulation (a dropped `ValueError` message and a dropped `MissingAdjoint`
substring) and are corrected here. The 16 group as:

| family | items | why it breaks |
|---|---|---|
| **G1** per-face wiring | 4 | `bc.apply(face_slot)` |
| **G3** transpose capability | 4 | `bc.apply_transpose(np.zeros_like(face_in))` |
| **G5** realizer-wiring probes | 8 (2 in `test_bound_compat.py`, 6 in `test_snmesh_realizer_wiring.py`) | assert the realized object IS a full-ordinate `PermutationOperator` / `TensorProductOperator`, and probe it with `np.ones(quad.N)` |

**G2** (the 2-D balance gate) and **G4** (the `_NoTransposeLaw` surrogate)
are NOT in the 16 — G2 needs no change at all, and G4 survives only because
its stub's `apply` is never called on the transpose path. G4 is still a
landmine (§D3) and belongs in the rewrite list.

> **C-1 verdict: the rewrite list is 16 breaking items across 4 files plus
> G2's prose and G4's stub — not "5 items across 2".** All are re-posed,
> none deleted.

---

## D2 — ⚠ THE CENTRAL QUESTION: can each gate still RED?

### D2.0 Method, and its positive control

The narrowing is applied in-process by `b3plugin.py`. Two fidelities:

* `ORPHEUS_B3=narrow` — narrowing at the **call site** only; `sn.bc[face]`
  keeps its full-face signature (so the old gates can still *run*, which is
  what lets their assertions be isolated).
* `ORPHEUS_B3=realizer` — **full fidelity**: `SNMesh.realize_boundary_law`
  emits a `Γ₊ → Γ₋` operator, so `sn.bc[face].apply` itself carries the
  narrowed signature. This is what B3 ships.

**[M] Positive control for every variant** (`b3/probe_variants.py` —
sha256 of `B.apply(ψ).boundary.values` and `B.apply_transpose(ψ)…`, four
fixtures):

| variant | `slab_vac_refl` FWD | `slab_refl_refl` FWD | `sphere_refl` FWD | `cyl_refl` FWD |
|---|---|---|---|---|
| **`none` (baseline)** | `ef8b015dd6855974` | `81dc734e71496d63` | `b4e9fa3dc07744fc` | `80d8b15e629155ee` |
| `narrow` | **identical** | **identical** | **identical** | **identical** |
| `realizer` | **identical** | **identical** | **identical** | **identical** |
| `realizer_lax` (reduced permutation, no guard) | **identical** | **identical** | **identical** | **identical** |
| `+N1` row-swap → outflow | `6bcb5442c56b6470` | `ad400b5061c6d965` | `c500f706e59e1f67` | `eee9e576f5226131` |
| `+N2` pass `Γ₊` through | `669234f8384fe64e` | `5b2e5f4018a96dd3` | `5dcf1c146deb9bb4` | `eaa9c207c407bc46` |
| `+N3` wrong half (`Γ₋` for `Γ₊`) | `49a11ff205fb79ae` | `3dc681398fdc2f03` | `85ba6fca4433c64f` | `181845adf282e9c2` |
| `+N4` full face (strict) | **RAISES** `ValueError` | RAISES | RAISES | RAISES |
| `+N4` full face (lax) | `49a11ff205fb79ae` | `c6cbfa23223162b1` | `85ba6fca4433c64f` | `e234881d17461484` |
| `+N5` tangent leak | identical | identical | identical | **`66b7617d9a8f0461`** (nnz 2→4) |
| `+N7` transpose row-swap (TRA column) | `49a11ff205fb79ae` | `3dc681398fdc2f03` | `85ba6fca4433c64f` | `181845adf282e9c2` |

Every mutation used below changed production numbers. **N5's earlier form was
VACUOUS** (a permutation's image on the tangent rows is structurally zero once
the domain is `Γ₊`) — caught by the control and re-spelled before any colour
was read. That is the `vv` Mode-8 METHOD WARNING firing in practice.

### D2.1 Per-gate verdict

| gate | survives the narrowing? | reddenable after? | by WHAT mutation |
|---|---|---|---|
| **G1 leg A** — `got[inflow] == bc.apply(face_slot)[inflow]` | **NO — it does not even RUN.** Its reference `bc.apply(psi.boundary.face_view(face))` feeds a full face to a `Γ₊`-domain operator | n/a — must be re-posed | — |
| **G1 leg B** — `assert not got[outflow].any()` | **YES, syntactically** — and **it CAN still red** | **YES** | **N2** (pass `Γ₊` through to the output outflow rows) and **N1** (write the image to the outflow rows). Both measured RED with the original message. |
| **G2** — 2-D balance | **YES, unchanged and green** (it never touches the law directly; it reads `B.apply(ψ)[inflow]`) | **YES** | **N1**, **N2**, **N3**, **N4(lax)** — all measured RED via the `keff ≠ k_inf` precondition |
| **G3** — `Bᵀ = lawᵀ ∘ P_inflow` | **NO — does not run** (same full-face reference, `:211`) | n/a — re-pose | — |
| **G4** — `_NoTransposeLaw` surrogate | runs, but its stub is a full-face identity — a latent surrogate landmine | yes, once the stub is narrowed | disable the `adjointable(law)` guard |
| **G5** — realizer-wiring probes | **NO — do not run** (`np.ones(quad.N)` probes) | n/a — re-pose | — |

### D2.2 The headline answer, with the measurement

**G1 leg B does NOT become a tautological guard.** Measured, `narrow+N2`:

```
[b3plugin] variant=narrow+N2  |B.apply|_max=1.845820701454e+00  nnz=6
FAILED tests/sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[slab_vacuum_reflective] - AssertionError: slab_vacuum_reflective face 'xmin': B emitted non-zero on t...
FAILED tests/sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[slab_reflective_reflective] - AssertionError: slab_reflective_reflective face 'xmin': B emitted non-zero ...
FAILED tests/sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[sphere_reflective] - AssertionError: sphere_reflective face 'xmax': B emitted non-zero on the ou...
FAILED tests/sn/operators/test_sn_boundary_operator.py::TestApply::test_apply_per_face_equals_legacy_bc_apply_on_inflow_row[cyl_reflective] - AssertionError: cyl_reflective face 'xmax': B emitted non-zero on the outfl...
8 failed, 9 passed, 1 warning in 1.37s
```

**But its TEETH CHANGE, and it acquires a hole.** After the narrowing:

* the bug it was written for — *the law's outflow image leaks into the
  output* — becomes **unspellable**: the law has no outflow image, its
  codomain is `Γ₋`. (Strictly better than a red; `vv` Mode-8's
  "unspellable ≠ red" note, lesson L24.)
* what remains reddenable is the **write-target** family: N1 (wrong target),
  N2 (extra write). Both measured.
* **it stays blind to the tangent rows** — the hole the review never
  recorded. Isolated measurement, `narrow+N5` (so the old gate can run):

```
=== old-gate isolation: narrow ===
8 passed, 20 deselected, 1 warning in 0.03s
=== old-gate isolation: narrow+N5 ===
FAILED ::test_rg2_codomain_is_exactly_the_inflow_trace[cyl_reflective]
1 failed, 7 passed, 20 deselected, 1 warning in 0.03s
```

The **old** `got[outflow]` leg is GREEN under a tangent leak; the **re-posed**
`got[outside Γ₋]` leg is RED. Positive control: `realizer+N5` moves the cyl
hash `80d8b15e629155ee → 66b7617d9a8f0461` and `nnz 2 → 4`.

> **So: keep the assertion, widen it.** `outflow` → `complement(inflow)`.
> The change costs nothing and closes a measured hole on the only gate
> fixture that has tangential ordinates.

### D2.3 The bug the narrowing INTRODUCES, and it is silent

**[M] `b3/probe_strict.py`** — a reduced `PermutationOperator` does **not**
validate its domain:

```
x_full  = [7.44   7.4635 5.1226 3.2864 1.4314 4.067  4.2678 1.3622]
x_half  = [7.44   5.1226 1.4314 4.2678]   (the honest Gamma_+)

PermutationOperator:
  honest Gamma_+ input  -> shape (4, 1) [4.2678 1.4314 5.1226 7.44  ]
  FULL-FACE input       -> NO RAISE, shape (4, 1) [3.2864 5.1226 7.4635 7.44  ]   equals-honest=False

perm & Identity:
  honest Gamma_+ input  -> shape (4, 1) [4.2678 1.4314 5.1226 7.44  ]
  FULL-FACE input       -> NO RAISE, shape (4, 1) [3.2864 5.1226 7.4635 7.44  ]   equals-honest=False
```

Fancy-indexing silently truncates. **Unless B3 ships an explicit domain
guard, "hand the law the full face" (N4) and "hand it the wrong half" (N3)
are both spellable and silent at the operator.** They are caught downstream
(§D3 RG-1, and G2's `keff` precondition), but the elegant fix is to make them
unrepresentable — see RG-3b.

---

## D3 — the re-posed contract, as assertions

All five are implemented and **measured** in
`$CLAUDE_JOB_DIR/tmp/b3/test_reposed.py`. They gate the **mechanism** (the
domain handed to the law, the codomain written), never a downstream aggregate.

### RG-1 — replaces G1 leg A: the law is fed `Γ₊`

```python
gamma_plus = psi.boundary.face_view(face)[trace.outflow_indices_for_face(face)]
np.testing.assert_array_equal(
    out.boundary.face_view(face)[inflow], bc.apply(gamma_plus),
    err_msg=f"{case_id} face {face!r}: B's inflow rows are not law(gamma_+) — "
            f"the law was handed the wrong domain.")
```

Reddens on: **N1** (row swap), **N3** (wrong half), **N4-lax** (full face),
and a face↔face swap (the asymmetric-BC slab param preserves that).

### RG-2 — replaces G1 leg B: the codomain IS `Γ₋`, stated as a WRITE-TARGET

```python
off_codomain = np.setdiff1d(np.arange(n_face), inflow)   # outflow  u  tangent
peak = float(np.abs(got[off_codomain]).max())
if peak != 0.0:
    pytest.fail(f"{case_id} face {face!r}: B emitted {peak:.3e} OUTSIDE its "
                f"codomain Gamma_- (rows {off_codomain.tolist()}) — the A_ss "
                f"block is not Gamma_+ -> Gamma_-.")
```

Reddens on: **N1**, **N2**, **N5** (tangent leak — the old leg's blind spot).
Uses `pytest.fail`, not a bare `assert`, so it is `-O`-proof in any module.

### RG-2b — the non-vacuity control RG-2 needs

RG-2 is satisfied by the **zero operator**. Pair it (never ship it alone):

```python
live = sum(np.count_nonzero(out.boundary.face_view(f)[inflow_f]) for f in faces)
if live == 0:
    pytest.fail(f"{case_id}: B emitted NOTHING on Gamma_- — RG-2 would be "
                f"vacuously satisfied by the zero operator.")
```

Reddens on **N1** (the image goes to the outflow rows, so `Γ₋` is empty).

### RG-3 — the SHAPE contract of the realized law (the honest new gate)

```python
probe = np.ones((outflow.size,) + tail)
image = sn.bc[face].apply(probe)
if image.shape != (inflow.size,) + tail:
    pytest.fail(f"realized law maps {probe.shape} -> {image.shape}; the "
                f"narrowed contract is Gamma_+ -> Gamma_-.")
```

This is the assertion C-1's *"state the narrowed contract"* actually means.
RED pre-B3, GREEN post-B3 — measured both ways.

### RG-3b — makes the D2.3 bug UNSPELLABLE

```python
with pytest.raises((ValueError, IndexError)):
    sn.bc[face].apply(np.ones(full_face_shape))
```

This is the **negative test for a guard B3 must add** — the reduced
permutation does not raise on its own (D2.3). Skips the face when
`n_outflow == n_face` (no distinguishable wrong shape), so it never
false-greens on a degenerate fixture.

### RG-4 — the transpose's codomain is `Γ₊`

```python
off = np.setdiff1d(np.arange(n_face), outflow)
if float(np.abs(got[off]).max()) != 0.0:
    pytest.fail("B^T emitted OUTSIDE its codomain Gamma_+ …")
```

Reddens on **N7** (transpose writes onto the inflow rows). Note the wide
transpose write B3 removes is **value-invisible** on every reachable law
(`+N6` measured sha-identical), so RG-4 is the only thing that pins it.

### RG-5 — NEW: the split partition on a **2-D** mesh (see §D5.4 for why)

```python
sn = _homogeneous_reflective_2d(nx=4, ny=4)          # 1-D is BLIND — measured
parts = SNBoundaryOperator(sn).split(SweepSchedule.gauss_seidel(sn))
total = parts.lower.apply(psi).boundary.values + parts.upper.apply(psi).boundary.values
np.testing.assert_array_equal(
    total, SNBoundaryOperator(sn).apply(psi).boundary.values,
    err_msg="B_lower + B_upper != B — the narrowed sel->position-within-Gamma_- "
            "remap does not partition the whole-trace reflection.")
```

Reddens on **N8**. The existing sibling
(`test_split_masked_halves_are_trace_only`, a **sphere**) stays green under
N8 — measured.

### G2 — the 2-D balance gate: keep the assertion, fix the vocabulary

G2 needs **no** assertion change (measured green under `narrow` /
`realizer`, and red under N1/N2/N3/N4). Its docstring/prose does — it
currently says *"``R·G·ψ.outflow`` is computed via the canonical
`SNBoundaryOperator`"* while the code passes the **whole field**. Post-B3
that sentence becomes true; note the change so a reader does not read the
old text as describing the old (wide) call.

### G4 — the duck-typed surrogate

`_NoTransposeLaw.apply(x) = x` must become
`apply(x) -> np.zeros((n_inflow,) + x.shape[1:])` (or a reduced identity when
`n_in == n_out`). **[R]** `tests/sn/operators/test_sn_boundary_operator.py:227`.
This is the FOURTH-search class from B2 — grep `SimpleNamespace(` and
`class _.*Law` near the boundary consumers before landing.

### V&V tagging for the re-posed set — and a clean bill on `catches`

**[G]** Neither gate file carries a `catches(...)` marker:

```
$ grep -rn "catches(" tests/sn/operators/test_sn_boundary_operator.py \
      tests/sn/operators/test_bc_extraction_2d.py \
      tests/geometry/test_bound_compat.py \
      tests/sn/operators/test_snmesh_realizer_wiring.py
(no output)
```

So B3 decays **no** ERR coverage claim — the L28 decayed-marker hazard does
not apply here. (Contrast B0.3, which had to relocate `catches("ERR-052")`.)

**[R]** `tests/sn/operators/test_sn_boundary_operator.py:44` —
`pytestmark = [pytest.mark.foundation]`, no `verifies(...)`. RG-1…RG-4 are
**software/structural invariants with no theory-page `:label:`**, so they
inherit `foundation` and MUST NOT carry `verifies(...)` (lesson L9 — stacking
them is silent level conflation).

**[R]** `tests/sn/operators/test_bc_extraction_2d.py:403` — G2's class carries
`@pytest.mark.l1` with **no** `verifies(...)`, so it shows as an orphan in the
audit. Out of B3's scope, but worth a one-line `verifies("…")` while the file
is open.

### The measured colour matrix for the re-posed set

```
################ ORPHEUS_B3=none ################
12 failed, 40 passed, 1 warning in 0.89s          <- RG-1/RG-3/RG-3b RED pre-B3, as designed
################ ORPHEUS_B3=realizer ################
8 failed, 44 passed, 1 warning in 0.89s           <- the 8 reds are the OLD G1+G3; ALL 24 re-posed items GREEN
################ ORPHEUS_B3=realizer+N1 ################
25 failed, 27 passed, 1 warning in 0.50s          <- RG-1, RG-2, RG-2b all RED
################ ORPHEUS_B3=realizer+N2 ################
17 failed, 35 passed, 1 warning in 1.01s          <- RG-2 RED
################ ORPHEUS_B3=realizer+N3 ################
15 failed, 37 passed, 1 warning in 0.50s          <- RG-1 RED
################ ORPHEUS_B3=realizer+N5 ################
10 failed, 42 passed, 1 warning in 0.95s          <- RG-2[cyl] RED  (the tangent hole)
################ ORPHEUS_B3=realizer+N7 ################
13 failed, 39 passed, 1 warning in 0.92s          <- RG-4 RED x4
################ ORPHEUS_B3=realizer_lax ################
4 failed, 20 passed, 1 warning in 0.04s           <- ONLY RG-3b RED: no domain guard
################ ORPHEUS_B3=realizer_lax+N4 ################
8 failed, 16 passed, 1 warning in 0.04s           <- RG-1 RED x4 + RG-3b RED x4
```

Every row above has its positive control in the §D2.0 hash table.

---

## D4 — the bit-identity gate

### D4.1 The narrowing IS value-neutral on every reachable SN configuration

**[M]** `b3/probe_narrow3.py` — real `SNMesh` fixtures, real `sn.bc[face]`,
real `angular_trace`; the question asked is *does `law.apply(full)[inflow]`
depend on anything outside `full[outflow]`?*

```
### slab_vacuum_reflective
  face=xmin  kind=vacuum       shape=(4, 1) in=2 out=2 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00
  face=xmax  kind=reflective   shape=(4, 1) in=2 out=2 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00

### slab_reflective_reflective
  face=xmin  kind=reflective   shape=(4, 1) in=2 out=2 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00
  face=xmax  kind=reflective   shape=(4, 1) in=2 out=2 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00

### sphere_reflective
  face=xmax  kind=reflective   shape=(4, 1) in=2 out=2 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00

### cyl_reflective
  face=xmax  kind=reflective   shape=(8, 1) in=2 out=2 tangent=4 NARROW-BIT-ID=True max|d|=0.000e+00

### cart2d_reflective
  face=xmin  kind=reflective   shape=(24, 2, 4) in=12 out=12 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00
  face=xmax  kind=reflective   shape=(24, 2, 4) in=12 out=12 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00
  face=ymin  kind=reflective   shape=(24, 2, 4) in=12 out=12 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00
  face=ymax  kind=reflective   shape=(24, 2, 4) in=12 out=12 tangent=0 NARROW-BIT-ID=True max|d|=0.000e+00
```

End-to-end, the whole tree under the honest narrowing:

```
2222 passed, 5 skipped, 2 xfailed, 28 warnings in 22.50s
```

(same 2222-item set the M8 sweep reported as `48 failed, 2174 passed`).

### D4.2 Design of the gate

* **Mesh × law × method rows.** The four `_CASES` fixtures
  (slab-asymmetric / slab-symmetric / sphere / **cyl with the product
  quadrature — mandatory, it is the only fixture with tangential
  ordinates**) plus the 2-D Cartesian `_homogeneous_reflective_2d(4,4)`
  (`ng=2`, the only multi-group + multi-spatial-axis face slot). Both
  `apply` and `apply_transpose`.
* **`np.array_equal`, not nulp.** The narrowing removes *rows* from a
  gather; it changes no reduction tree. Measured: sha256 identical on all
  four fixtures for BOTH `_NarrowedLaw` (a wrapper) and `_LaxNarrowedLaw`
  (an **independently written** reduced-gather implementation). A
  `nulp` tolerance here would hide exactly the class of bug the gate exists
  for. Per `vv` §bit-identity: this is criterion-3 *zero* drift, so the
  contract stays `array_equal`.
* **Capture the baseline BEFORE the carve** — the campaign's
  `--capture-baseline` conftest flag, or a committed sha256 of
  `B.apply(ψ).boundary.values` / `B.apply_transpose(ψ)…` per fixture (the
  eight hashes in §D2.0 are already the pre-carve values at `acb46245`).
  Old-vs-new is **necessary, not sufficient** (lesson L2) — pair it with the
  existing closed-form anchors that already sit downstream: G2's
  `keff == k_inf == 1.875` and
  `test_streaming_collision_operator.py::…::test_fixed_source_homogeneous_reflective_recovers_q_over_sigma_composite`
  (`Q/Σ_t`). Both are in the 43-item consumer set that must stay green.

### D4.3 The seven laws — mesh-level for two, LAW-level for five

**[M]** SN admits `{reflective, vacuum}` only, so a mesh-level gate cannot
reach the other five. And they do **not** all narrow neutrally
(`b3/probe_narrow.py`, operator-level, `SNMethodSpace` + realizer):

| law | narrows bit-identically? | why |
|---|---|---|
| `reflective` (α=1, α=0.7) | **YES** | `perm[inflow] ⊆ outflow` — measured `subset-of-outflow=True` on gl4 / lebedev5 / LS(4) |
| `vacuum` | **YES** | the inflow rows of the image are identically zero |
| `white` (α=1, α=0.4) | **YES** at the correctly-oriented face | `cos_w` is supported on the **outgoing** hemisphere only |
| `prescribed_inflow` | **YES** | a constant source, input-independent |
| **`albedo(α∉{0,1})`** | **NO** | **[R]** `sn/boundary/realizer.py:312-326` realizes it as `α·(IdentityOperator() & IdentityOperator())` — a **full-face** identity, so its inflow-row output is `α·γ₋ψ`. Narrowing changes it (measured `max|a-b|` 0.44–0.67 on every quadrature). |
| **`periodic`** | **NO** | same shape — the realized `PeriodicWrapOperator` is identity-like (issue **#183** says its `apply` is `x.copy()`) |
| `zero_flux` | n/a | SN refuses it (`BoundaryError`) |

> **This is a finding, not just a gating note.** The SN-realized albedo law
> re-emits the **incoming** flux scaled by α, not the outgoing one. It is
> unreachable today (`{reflective, vacuum}`), so nothing is wrong in
> production — but **#189 / B6 make it reachable**, and B3 is the phase that
> decides its `Γ₊ → Γ₋` action. Do not let the narrowing silently "fix" it:
> land it as a deliberate, stated change with its own value gate, or B6
> inherits a law whose pre-B3 and post-B3 meanings differ.

**Recommended gate shape:**

* **Mesh-level `array_equal`** for `{reflective, vacuum}` × 5 fixtures ×
  {`apply`, `apply_transpose`} — the bit-identity claim proper.
* **Law-level** for the other five: assert the realized narrowed operator's
  `(domain, codomain)` shapes (RG-3) and, for `white` / `prescribed_inflow`,
  `array_equal` against `full_law.apply(embed(γ₊))[inflow]` computed from the
  **pre-B3 realizer** — a law-level bit-identity that needs no mesh. For
  `albedo` and `periodic` a bit-identity claim is **false**; ship them as
  `xfail(strict=True, reason="B3: albedo/periodic realize as a full-face
  identity; their Γ₊→Γ₋ action is decided in B4/B6 — #183, #189")`, which is
  the campaign's own todo-list technique and makes the decision impossible to
  forget.

---

## D5 — the mutation sweep

### D5.1 Where the set lives

`$CLAUDE_JOB_DIR/tmp/` (this job, `c30e4f25`) — the B0.3 notes
(`scratch/b0_3_coverage_repairs.md:29-45`) name it and it is still on disk:

| file | mutations | env var |
|---|---|---|
| `mutplugin.py` | 12 leaf-action | `ORPHEUS_MUT` |
| `mutplugin2.py` | 14 guard-disabling | `ORPHEUS_GUARD` |
| `mutplugin3.py` | 5 adjoint | `ORPHEUS_ADJ` |
| `mutplugin4.py` | the ERR-052 probe | `ORPHEUS_ERR052` |
| `run_mutations.sh` / `run_guards.sh` / `run_adj.sh` | runners (17-file target list) | — |

**It is NOT in the repo.** A campaign whose stated gate is *"the mutation
sweep still catches 30/31"* (plan §B3) currently depends on a job-scratch
directory that evaporates. **B3 should promote the harness to a tracked
location** (`tests/_harness/mutations/` or `scratch/`), or the gate cannot be
re-run by the next session.

### D5.2 Which mutations B3 breaks or makes vacuous

| # | mutation | fate under B3 | action |
|---|---|---|---|
| **M8** | `reflect_nosel` — drop the `sel` restriction (`out[...] = law.apply(face_in)`) | **UNSPELLABLE.** There is no full-face image to write; `law.apply(γ₊)` has shape `(n_inflow, …)` and `out.face_view(face)[...]` has `(n_face, …)`. The monkeypatch will `ValueError`, i.e. the harness aborts rather than measuring | **replace with N2** (write `γ₊` through onto the outflow rows) — the B3 analogue, measured RED on G1 leg B and RG-2 |
| **M6** | `vacuum_complement` — rebind `IncomingOrdinateMaskTensor` in the realizer | **AT RISK of going vacuous.** If B3 realizes vacuum as a rectangular zero map `Γ₊→Γ₋`, the mask tensor is no longer on the path (the review §5 predicted exactly this collapse of "two spellings of `P_in` + one of `I − P_in`" into one projector). The measured `_LaxNarrowedLaw` vacuum arm never touches it | **re-point** at whatever B3 emits for vacuum, and re-run the bite check |
| **M7** | `psrc_nomask` — drop `IncomingSourceOperator`'s `Γ₋` mask | **AT RISK.** After narrowing, the source law's codomain IS `Γ₋`, so an explicit mask is structurally redundant → the mutation becomes a no-op | **re-point** or retire with a stated reason |
| **M1 / M2** | reflective perm → identity / rolled | **survive**, but now act on the **reduced** permutation. `np.roll(perm,1)` on a reduced index set is a different (still valid) index drift | **re-verify the bite**; the hand-computed catchers in `test_sn_boundary_realizer.py` also feed `np.ones(quad.N)` and are in the §D1.3 rewrite list |
| **M3 / M4 / M5** | white weights | survive at the operator level (white is SN-unreachable either way) | re-verify if `AngularAverageOperator.from_quadrature` gains a `Γ₊` signature |
| **M9** | SN `α → 1−α` at the realizer seam | **survives unchanged** | — |
| **M10 / M11 / M12** | diffusion | **untouched** — diffusion is already narrow | — |
| guard set (14) | law invariants + realizer refusals | **untouched** — they fire before realization | — |
| adjoint set (5) | `.H` / reciprocity | **survive by value** (the transpose is sha-identical) but any that build a full-face probe are in the rewrite list | re-verify |

### D5.3 The NEW mutations B3 owes

Your instinct is right, and it is now measured. Six rows, all with a
positive control in §D2.0:

| id | mutation | catcher | measured |
|---|---|---|---|
| **N4** | **hand the law the FULL face** instead of `γ₊` | **RG-3b** (the domain guard's negative test) if B3 ships the guard; otherwise **RG-1** ×4 and **G2**'s `keff` precondition | `realizer+N4` RAISES with the guard; `realizer_lax+N4` gives `keff = 1.767 ≠ 1.875` and RG-1 RED ×4 |
| **N3** | hand the law the **wrong half** (`γ₋` for `γ₊`) — spellable and silent *because `n_inflow == n_outflow` on every reachable fixture* | **RG-1**, **G2** | `realizer+N3` → 15 failed incl. RG-1 ×4 |
| **N1** | write the image to the **outflow** rows | **RG-1**, **RG-2**, **RG-2b**, G1 leg B, G2 | `realizer+N1` → 25 failed |
| **N2** | pass `γ₊` **through** to the output outflow rows (the M8 replacement) | **RG-2**, G1 leg B, G2 | `realizer+N2` → 17 failed |
| **N5** | leak onto the **tangent** rows | **RG-2 only** — the old gate is blind | `narrow+N5`: old gate 8 passed, RG-2[cyl] RED |
| **N7** | transpose writes onto the **inflow** rows | **RG-4 only** (plus the dense-`Bᵀ` oracle in `test_psi_half_coupling.py`) | `realizer+N7` → 13 failed incl. RG-4 ×4 |
| **N8** | `pos = arange(sel.size)` instead of `searchsorted(inflow, sel)` on the split path | **nothing in the boundary suite** — only `tests/sn/solve/test_gauss_seidel_reification.py` + three 2-D solve tests. **Needs a new 2-D split gate** (§D5.4) | 1-D sha-identical (vacuous); 2-D lower/upper both move, `lower+upper ≠ whole` |

### D5.4 ⭐ The mutation the plan could not have predicted — N8, and its 1-D blind spot

The narrowing introduces a **brand-new piece of index arithmetic** that does
not exist today. Currently the split path writes `full[sel]` — `sel` indexes
the *face*, and `full` is a *face-shaped* image, so no remap is needed. After
narrowing, the law's image is indexed by **position within `Γ₋`**, so the
Gauss-Seidel row subset must be remapped:

```python
pos = np.searchsorted(inflow, sel)          # NEW — a Mode-5 hazard that did not exist
out_boundary.face_view(face)[sel] = image[pos]
```

**N8** = the plausible transcription `pos = np.arange(sel.size)` (i.e. assume
`sel` is a prefix of `inflow`).

**[M] N8 is VACUOUS in 1-D and BITES in 2-D** — because in 1-D the schedule
hands each face *entirely* to one half:

```
=== ORPHEUS_B3=none ===                        === ORPHEUS_B3=realizer+N8 ===
  SLB lower: sha=c4852f936cb2cc11                SLB lower: sha=c4852f936cb2cc11   <- unchanged
  SLB upper: sha=9f3e6044b46abe3b                SLB upper: sha=9f3e6044b46abe3b   <- unchanged
  CYL upper: sha=f2f5fc87969ca9cc                CYL upper: sha=f2f5fc87969ca9cc   <- unchanged

  2D lower: sha=ca03a7cf526cdd9c                 2D lower: sha=cc5e5d07550257a8    <- CHANGED
  2D upper: sha=41c99ea99a6185f8                 2D upper: sha=ce901612dbc0e885    <- CHANGED
  2D whole: sha=bdc5be4b606bfd46                 lower+upper==whole: False
```

The activating structure, printed from the live schedule:

```
   face ymin: inflow=[2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23]
      lower rows=[6, 7, 14, 15, 22, 23]        <- NOT a prefix of inflow
      upper rows=[2, 3, 10, 11, 18, 19]
   face xmin: inflow=[4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23]
      lower rows=[4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23]   <- the WHOLE face (1-D-like)
      upper rows=[]
```

**[M] And the boundary suite's own split gate is BLIND to it** — it is built
on a **sphere** (`_sphere(bc="reflective")`,
`tests/sn/operators/test_psi_half_coupling.py:795`), a 1-D mesh. Under
`realizer+N8`:

```
FAILED tests/sn/solve/test_gauss_seidel_reification.py::test_w2_split_exactness
FAILED tests/sn/solve/test_gauss_seidel_reification.py::test_w2_fixed_point_equivalence_diagonal_cubature
FAILED tests/sn/solve/test_2d_anisotropic_windowing.py::test_2d_p1_aniso_moment_path_carries_signal_and_si_krylov_agree
FAILED tests/sn/solve/test_affine_carve_bit_identity.py::test_converged_flux_bit_identical_after_affine_carve[si_2d_p1_aniso_het]
FAILED tests/sn/solve/test_fixed_source_2d_equivalence.py::test_2d_heterogeneous_si_krylov_equivalence
21 failed, 2324 passed, 4 skipped, 3 deselected, 2 xfailed, 14 warnings in 174.88s (0:02:54)
```

`test_split_masked_halves_are_trace_only` is **not** in that list — it stayed
GREEN. The catchers all live in `tests/sn/solve/`, i.e. **outside** the
boundary suite, and they are end-to-end solves rather than mechanism gates.

> **B3 owes a 2-D split gate in the boundary suite** — the `§0.6`
> config-blindness rule in its purest form: the convenient (1-D) fixture nulls
> the exact term the carve introduces. Concretely: run
> `test_split_masked_halves_are_trace_only` on
> `_homogeneous_reflective_2d(4,4)` **in addition to** the sphere, asserting
> `lower + upper == whole` bit-identically. Measured RED under N8, GREEN under
> `realizer`.

Plus one **negative** row that must stay GREEN and is the honest control:

| id | mutation | expectation | measured |
|---|---|---|---|
| **N6** | the transpose keeps the OLD wide (full-face) write target | **GREEN — value-invisible** on every reachable law | `narrow+N6` sha-identical; 17 passed |

N6 is why RG-4 is needed: without it, "the transpose's codomain narrowed"
is an unverifiable claim.

---

## The re-posed gate set, as a checklist

| id | replaces / adds | assertion | reddens on | fixture requirement |
|---|---|---|---|---|
| **RG-1** | G1 leg A | `got[inflow] == bc.apply(face[outflow])` | N1, N3, N4-lax, face↔face swap | keep the asymmetric-BC slab param |
| **RG-2** | G1 leg B | `got[complement(inflow)] == 0`, via `pytest.fail` | N1, N2, **N5** | **`cyl_reflective` is mandatory** — the only fixture with tangential ordinates |
| **RG-2b** | new (non-vacuity control for RG-2) | `Γ₋` rows are not all zero | N1 | any |
| **RG-3** | new | realized law maps `(n_out, …) → (n_in, …)` | pre-B3 state | all four |
| **RG-3b** | new (makes N4 unspellable) | `pytest.raises` on a full-face input | absence of the domain guard | skip when `n_out == n_face` |
| **RG-4** | new | `Bᵀ`'s image is supported on `Γ₊` alone | N7 | all four |
| **RG-5** | new (§D5.4) | `B_lower + B_upper == B` on a **2-D** mesh, `array_equal` | **N8** | `_homogeneous_reflective_2d(4,4)` — 1-D is BLIND |
| **G2** | unchanged assertion | `defect < 10 × solver_tol` | N1, N2, N3, N4 | fix the docstring's "R·G·ψ.outflow" wording |
| **G3** | rewrite reference | `Bᵀ = lawᵀ ∘ P_Γ₋` with a `Γ₋`-shaped `masked` | disabling the `adjointable(law)` guard | all four |
| **G4** | narrow the stub | `_NoTransposeLaw.apply(x) -> zeros((n_in,) + x.shape[1:])` | same | — |
| **G5** | rewrite probes | probe with `np.ones((n_out,) + tail)`; assert the reduced structure | M1/M2 (re-verified on the reduced perm) | 8 items |
| **bit-id** | new | `array_equal` on `B.apply` / `B.apply_transpose`, 5 fixtures, pre-vs-post | any value change | + law-level rows for the 5 unreachable laws |

---

## Summary of findings the plan should absorb

1. **C-1's rewrite list is 16 breaking items across 4 files (+ G2's prose and
   G4's stub), not 5 items across 2** — the three extra families are the
   transpose-capability tests, the `_NoTransposeLaw` duck-typed surrogate, and
   the eight realizer-wiring probes. (§D1.3)
2. **G1 leg B does not go tautological** — it stays reddenable by N1/N2. But
   its original tooth (the law's outflow image leaking) becomes
   **unspellable**, and it is **blind to the tangent rows**. Widen
   `outflow` → `complement(inflow)`. (§D2.2)
3. **The face slot is a THREE-way split** — `inflow ⊔ outflow ⊔ tangent`,
   and the `cyl_reflective` gate fixture has 4 tangential ordinates out of 8.
   Every "the rows are zero" assertion in this subsystem should be audited
   against that. (§0, §D2.2)
4. **The narrowing is bit-identical on every reachable SN configuration** —
   sha256-identical forward AND transpose on five fixtures, under two
   independent implementations of the narrowed law; `2222 passed / 0 failed`
   whole-tree. The gate is `array_equal`, never nulp. (§D4)
5. **A reduced `PermutationOperator` does not validate its domain** — so
   B3's own headline bug is silent unless the carve adds a guard. RG-3b is
   its negative test. (§D2.3)
6. **`albedo` and `periodic` do NOT narrow neutrally** — they realize as
   full-face identities, so their `Γ₊→Γ₋` action is a decision B3 makes, not
   a refactor it inherits. Unreachable today; **#189 / B6 make them
   reachable**. Ship them as strict xfails so the decision cannot be
   forgotten. (§D4.3)
7. **The 31-mutation harness is not in the repo.** The B3 gate "the sweep
   still catches 31/31" is unrunnable by the next session until it is
   promoted to a tracked path. (§D5.1)
8. **M8 becomes unspellable; M6 and M7 are at risk of going vacuous.**
   Re-verify every leaf mutation's bite AFTER the carve — a mutation that
   stops biting silently deletes coverage. (§D5.2)
9. **B3 introduces a NEW index remap (`sel` → position within `Γ₋`) whose
   only activating regime is 2-D**, and the boundary suite's split gate is
   1-D (spherical) — measured blind. The catchers today are end-to-end 2-D
   solves in `tests/sn/solve/`. Add RG-5, a 2-D mechanism gate. (§D5.4)
