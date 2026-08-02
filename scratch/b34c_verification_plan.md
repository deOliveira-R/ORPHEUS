# B3.4c — verification plan (periodic narrows to the partner face)

Branch `refactor/operator-strategy-layers`, HEAD `cea73f8d`. Written by
test-architect 2026-08-01. Design of record:
`.claude/plans/b3_domain_narrowing_crosswalk.md` §11.2 (the periodic ruling),
§11.3 (SCC), §11.5 (doc falsehoods), §16.4 (the periodic gate already written
and red), §16.5 (the recurring index-order hazard). Blast radius already mapped
by explorer at `scratch/b34c_periodic_blast_radius.md`.

Everything marked `[M]` was MEASURED in this session; verbatim stdout is quoted
at the point of claim.

> ⚠ **L-007 intra-session drift is ACTIVE.** The working tree moved WHILE this
> plan was being written. As of the final read, steps **1, 2, 3, 5 and 6 of the
> plan of record are ALREADY IMPLEMENTED, uncommitted**, and step 4 is not. See
> §0.0 for the exact state; re-derive before acting.

---

## 0.0 `[M]` What has already landed (uncommitted working tree)

```
$ git status --porcelain -- orpheus/ tests/
 M orpheus/geometry/boundary/_factors.py
 M orpheus/numerics/face_layout.py
 M orpheus/numerics/spaces/angular_trace_space.py
 M orpheus/sn/boundary/angular.py
 M orpheus/sn/boundary/realizer.py
 M orpheus/sn/loss_representation/sweep_schedule.py
 M orpheus/transport/method.py
 M tests/numerics/test_angular_trace_space.py
```

| plan step | state | evidence |
|---|---|---|
| 1 — mint the bijection in `face_layout.py`, repoint 5 sites | **DONE** | `face_name` / `face_normal` / `FACE_NAMES` + `_FACE_SUFFIX_SIGN` / `_FACE_SIGN_SUFFIX` (+156 lines); all five sites read them |
| 2 — `opposite_face` | **DONE**, named `face_opposite` | `face_layout.py` |
| 3 — `BoundaryGeometryMap.domain_face(face)` + `SpatialWrap` axis guard | **DONE** | `_factors.py` (+128); `IdentityMap`/`SpecularMirror` return `face`, `SpatialWrap` returns `face_opposite(face)` and RAISES on an axis mismatch |
| 4 — `_reflect_trace` restricts at the SOURCE face | **NOT DONE** | `orpheus/sn/operators/boundary.py` is clean; the realizer's new comment already names the unbuilt consumer `SNBoundaryOperator._face_domains` |
| 5 — assert `Γ₊(partner) ≡ Γ₋(face)` at realization | **DONE** | `realizer.py::_assert_wrap_identification` |
| 6 — `SpatialWrap.is_adjointable` `False → True` | **DONE** | `_factors.py` |
| 7 — five doc falsehoods | **PARTIAL** — `_factors.py` and `angular.py` repaired; `numerics/operator.py:2616-2641` (`PeriodicWrapOperator`, falsehoods 1 + 3) and `geometry/boundary/periodic.py:35-45` (falsehood 2) still stand verbatim | Read |

`[M]` The live tree's damage, measured on the eleven affected files:

```
$ .venv/bin/python -O -m pytest tests/geometry tests/numerics/test_periodic_wrap_operator.py \
    tests/numerics/test_operator_capability_predicates.py tests/numerics/test_face_layout.py \
    tests/numerics/test_angular_trace_space.py tests/sn/operators/test_b3_domain_narrowing.py \
    tests/sn/operators/test_sn_boundary_realizer.py tests/sn/operators/test_operator_block_role.py \
    tests/sn/operators/test_boundary_conditions.py tests/sn/operators/test_snmesh_realizer_wiring.py \
    tests/diffusion/test_boundary_realizer.py -q
FAILED tests/geometry/test_boundary.py::test_periodic_bc_returns_input_unchanged - orpheus.geometry.boundary._errors.BoundaryError: SNBoundaryRealizer cannot ...
FAILED tests/sn/operators/test_sn_boundary_realizer.py::TestRealizePeriodic::test_periodic_returns_tensor_product_and_passes_through - orpheus.geometry.boundary._errors.BoundaryError: SNBoundaryRealizer cannot ...
2 failed, 782 passed, 4 skipped, 5 xfailed, 5 warnings in 5.87s
```

Both reds are the same cause, and it is the intended one:

```
E  orpheus.geometry.boundary._errors.BoundaryError: SNBoundaryRealizer cannot realize
   'periodic' without outflow_indices: since B3.2 a boundary law is typed Γ₊ → Γ₋, and its
   DOMAIN is the outflow half-trace of a particular face — which a quadrature alone cannot
   name. Construct via SNMethodSpace.for_face(quadrature=..., face=..., trace=...) ...
orpheus/sn/boundary/realizer.py:224: BoundaryError
```

---

## 0.1 Executive summary — the five findings that change the plan

1. **`periodic` is NOT in `SNMesh.BOUNDARY_OPERATOR_REGISTRY`**
   (`orpheus/sn/mesh/augmented_mesh.py:175-178` = `{"vacuum", "reflective"}`), so
   `BC("periodic")` RAISES at mesh construction and **no tag-constructed
   `SNBoundaryOperator` can ever hold a periodic face**. Step 4 therefore lands
   as *unreachable-from-a-tag* production code unless #189 lands with it. The
   end-to-end gate is still buildable and I built it: `[M]`
   `sn.realize_boundary_law(PeriodicBoundary(axis="x"), face)` bypasses the tag
   registry and works today. **§C4 uses that injection; the plan owes an
   explicit decision on whether B3.4c also admits the tag.**
2. **The already-red snapshot gate cannot be made honest by the documented
   edit.** §16.4 prescribes "hand the operator the PARTNER's `Γ₊`". That edit
   makes `test_delivers_the_partner_faces_outflow` a **character-for-character
   duplicate** of the live sibling `test_the_identity_body_is_correct_on_the_partner_half_trace`
   — a row that is green TODAY and green after B3.4c. Deleting the xfail would
   signal "B3.4c landed" while asserting nothing new: vv Mode-8 class 4, a
   marker whose flip is a no-op. §B gives the honest edit and the exact lie to
   avoid.
3. **`SpatialWrap.is_adjointable: False → True` changes NO observable in the SN
   stack.** `[M]` a realized periodic face reports `is_adjointable == True`
   TODAY, because the composite predicate reads the realized
   `PeriodicWrapOperator` (identity body ⇒ `True`) and never the law factor —
   which is doc falsehood #4 stated as a measurement. So the step-6 gate must
   be on the **factor**, plus on a NET-NEW transpose that actually crosses
   faces, or it is designed-green. §C6.
4. **The discriminating quadratures are `gauss_legendre(*)` and `lebedev(17)`.**
   `[M]` on `product(2,4)`, `level_symmetric(4)` and `level_symmetric(6)` the
   mirror image of `Γ₋` coincides with `sorted(Γ₊)` (§16.5's hazard, measured
   below), so an identity-vs-permutation confusion is invisible there.
5. **Two MISATTRIBUTED strict xfails are already in the tree and name
   periodic.** `[M]` `tests/geometry/test_boundary.py::test_wave0_sum_of_realized_bcs_acts_as_weighted_sum`
   and `::test_operator_sum_of_bcs_acts_as_weighted_sum` xfail on a
   `BoundaryError` out of `_outflow_restriction` — a faceless
   `SNMethodSpace.minimal` in the test's own helper — **not** on the documented
   "mismatched factor shapes". Their reason text says "unstateable until B3.4
   narrows white / albedo / periodic", which is false on all three counts now.
   They will NOT flip when B3.4c lands. §A row A-22.
6. **The existing block-diagonality gates will NOT red — they will DECAY.**
   `[M]` all three (`test_sn_boundary_operator.py:332`, `:656-700` ×2) are posed
   on `(BC.reflective, BC.reflective)`, so B3.4c cannot touch them; what becomes
   false is their stated JUSTIFICATION ("``B`` is block-diagonal over faces, so
   the subset MUST be the EXACT restriction"). A green gate whose reason has
   decayed is vv Mode-8 class 7 and nothing notices. §C7 carries the rescoping
   list and the positive cross-face gate that replaces the claim.
7. **The newly-minted bijection is currently UNGATED.** `[M]` `face_name` (the
   function), `face_opposite` and `FACE_NAMES` have zero test coverage, and the
   docstring `>>>` examples do not run (`pyproject.toml:80` sets no
   `--doctest-modules`). §C0-C2.

---

## 0.2 ⚠ RECONCILIATION — step 4 LANDED mid-session; the gates were dry-run against it

Between the §0.0 snapshot and the end of this session the main agent implemented
**step 4** (`SNBoundaryOperator._face_domains` + `_reflect_trace` reading the
source face, `orpheus/sn/operators/boundary.py` +89/−8) and repaired the two
`test_boundary.py` mixed-law xfails (§A-22). `[M]` the tree at final read:

```
$ git status --porcelain -- orpheus/ tests/
 M orpheus/geometry/boundary/_factors.py       M orpheus/sn/operators/boundary.py
 M orpheus/numerics/face_layout.py             M orpheus/sn/solver.py
 M orpheus/numerics/spaces/angular_trace_space.py   M orpheus/transport/method.py
 M orpheus/sn/boundary/angular.py              M tests/geometry/test_boundary.py
 M orpheus/sn/boundary/realizer.py             M tests/numerics/test_angular_trace_space.py
 M orpheus/sn/loss_representation/__init__.py  M tests/sn/operators/test_sn_boundary_realizer.py
 M orpheus/sn/loss_representation/sweep_schedule.py
```

The landed `_face_domains` also carries a **permutation certification** over the
boundary faces (`boundary.py:334-347`) — stronger than anything §C proposed, and
it subsumes two of my negatives (half-declared pair, absent partner) at the
COMPOSITE rather than at realization. §C8 is re-posed accordingly below.

### 0.2.1 `[M]` Dry-run of every §C gate against the landed tree

The plan's gates were transcribed verbatim into
`scratch/b34c_proposed_gates_dryrun.py` and run. FIRST pass:

```
$ .venv/bin/python -O -m pytest scratch/b34c_proposed_gates_dryrun.py -q -p no:cacheprovider
2 failed, 93 passed, 2 warnings in 1.00s
FAILED …::test_B_is_block_structured_not_block_diagonal - AssertionError:
FAILED …::test_a_wrap_whose_partner_face_is_absent_is_REFUSED - Failed: DID NOT RAISE
```

After the two corrections below:

```
$ .venv/bin/python -O -m pytest scratch/b34c_proposed_gates_dryrun.py -q -p no:cacheprovider
95 passed, 1 warning in 0.74s
```

Every cross-face row now PASSES — step 4 works:

```
TestPeriodicIsACrossFaceCoupling::test_the_two_faces_carry_different_data       PASSED
TestPeriodicIsACrossFaceCoupling::test_the_inflow_is_the_PARTNERS_outflow[xmin-xmax] PASSED
TestPeriodicIsACrossFaceCoupling::test_the_inflow_is_the_PARTNERS_outflow[xmax-xmin] PASSED
TestPeriodicIsACrossFaceCoupling::test_it_is_NOT_the_faces_own_outflow[xmin]    PASSED
TestPeriodicIsACrossFaceCoupling::test_it_is_NOT_the_faces_own_outflow[xmax]    PASSED
TestPeriodicIsACrossFaceCoupling::test_nothing_lands_off_gamma_minus[xmin]      PASSED
TestPeriodicIsACrossFaceCoupling::test_nothing_lands_off_gamma_minus[xmax]      PASSED
test_euclidean_reciprocity_on_a_periodic_slab[4]                               PASSED
test_euclidean_reciprocity_on_a_periodic_slab[8]                               PASSED
test_the_transpose_scatters_over_the_PARTNERS_outflow                          PASSED
test_the_composition_supplies_the_partner_half_trace                           PASSED
```

**Two corrections the dry-run forced:**

1. **`test_B_is_block_structured_not_block_diagonal` must use a ZERO base.**
   The differencing form `B(bumped) − B(base)` is NOT bit-exact — `[M]` 2/4
   entries land at `0.9999999999999998` / `1.0000000000000002` (1 ULP, the
   catastrophic-cancellation of `(x+1) − x`). Exploiting linearity instead is
   bit-exact:

   ```
   ZERO-BASE form: xmin Gamma- = [1.]  bit-exact 1.0? True
     xmin off-Gamma- all zero? True
     xmax face entirely zero? True   <- xmax's Gamma- reads xmin's Gamma+, which is 0 here
   ```

   so the gate is `base = zeros`, `bumped = zeros with Γ₊(xmax) = 1.0`, then
   `assert_array_equal(B(bumped).face_view("xmin")[Γ₋], 1.0)` — and the third
   line above is a free extra leg (the OTHER face stays zero, which pins the
   direction).

2. **§C8 re-poses.** `[M]` realization is still silent on a sphere; the refusal
   now lives at `SNBoundaryOperator.__init__`-time via the permutation
   certification. So the gate becomes:

   ```python
   def test_a_wrap_whose_partner_face_is_absent_is_REFUSED_at_the_COMPOSITE():
       r"""A one-faced mesh is not a torus — and the refusal's HOME moved.

       `[M]` `sn.realize_boundary_law(PeriodicBoundary(), "xmax")` on a sphere
       still succeeds: the realizer's identification guard only compares the
       partner's Γ₊ against this face's Γ₋, which it can compute from the face
       NAME whether or not the mesh carries that face. The composite is where
       the mesh's face inventory is known, and `_face_domains`'s permutation
       certification is what catches it. Sited here so the diagnosis names the
       topology instead of surfacing as a KeyError in `_reflect_trace`.
       """
       sn = _sphere()                      # faces == ("xmax",)
       sn.bc["xmax"] = sn.realize_boundary_law(PeriodicBoundary(axis="x"), "xmax")
       with pytest.raises(ValueError, match="domain map"):
           SNBoundaryOperator(sn)._face_domains
   ```

   ⚠ **Open design question for the plan of record:** should the *realizer*
   also refuse, or is the composite the right home? Composite-only means a
   realized wrap on a one-faced mesh is a legal object that cannot be composed
   — defensible (the law is fine, the CONFIGURATION is not), but it must be
   stated, because "realize succeeded" reads as "this is installable".

### 0.2.2 `[M]` Mutation battery — every gate's teeth, measured

In-process `monkeypatch` autouse plugins; **no `git checkout`** (the tree carries
uncommitted work). Baseline `NONE` = 3 failures: the two above plus one harness
artifact (a relative snapshot path resolved against the runner's cwd).

```
NONE                                 95 passed, 1 warning in 0.84s
M1_suffix_sign_swapped               12 failed, 83 passed
M3_face_opposite_is_identity         18 failed, 77 passed
M4_domain_face_is_self               25 failed, 70 passed
M8_reflect_reads_own_face             7 failed, 88 passed
M12_wrap_not_adjointable              2 failed, 93 passed
M13_wrap_body_scales                  6 failed, 89 passed
```

`[M]` per-mutation, the EXACT rows that reddened (`-rf`, names elided of the path):

```
### M8_reflect_reads_own_face   (revert step 4: _face_domains -> {f: f})
    TestPeriodicIsACrossFaceCoupling::test_the_inflow_is_the_PARTNERS_outflow[xmin-xmax]
    TestPeriodicIsACrossFaceCoupling::test_the_inflow_is_the_PARTNERS_outflow[xmax-xmin]
    TestPeriodicIsACrossFaceCoupling::test_it_is_NOT_the_faces_own_outflow[xmin]
    TestPeriodicIsACrossFaceCoupling::test_it_is_NOT_the_faces_own_outflow[xmax]
    test_the_transpose_scatters_over_the_PARTNERS_outflow
    test_B_is_block_structured_not_block_diagonal
    test_a_wrap_whose_partner_face_is_absent_is_REFUSED_at_the_COMPOSITE

### M12_wrap_not_adjointable
    test_the_wrap_declares_an_honest_transpose
    test_every_deck_transformation_is_adjointable

### M1_suffix_sign_swapped
    test_the_bijection_reproduces_the_retired_transcriptions[0--1] … [2-1]   (6)
    test_facelabel_render_agrees_with_the_bijection[0-min--1] … [2-max-1]    (6)

### M3_face_opposite_is_identity
    test_the_wrap_consumes_the_OPPOSITE_face[x] [y] [z]
    test_the_identification_compares_SETS_not_sizes
    test_the_honest_space_realizes
    (…all 7 TestPeriodicIsACrossFaceCoupling / reciprocity / transpose /
      block-structure / composite-refusal / composition rows)
```

**Three predicted blindnesses CONFIRMED by measurement, not asserted:**

* **M8 does NOT red `test_euclidean_reciprocity_on_a_periodic_slab`** — the
  reciprocity gate is blind to the B3.4c defect exactly as §C5/§D claim.
* **M3 does NOT red `test_opposite_is_an_involution`** — the identity IS an
  involution; `test_opposite_has_no_fixed_point` is the catcher (§E.1).
* **M1 does NOT red either round-trip row** — a bijection's round-trip is
  invariant under relabelling both directions; the `legacy` transcription leg
  and the `FaceLabel` agreement leg are the catchers (§E, M1).

**M12 reds exactly two rows and nothing else** — the honest scope of §C6: the
factor flip has no downstream observable.

**M6 (identification by SIZE) is proven separately** by probe F (§E.2), not by a
mutant: `[M]` the real guard raises on a same-sized-wrong-set space and a size
comparison would pass. A substitute-function mutant also broke the control leg,
so the direct measurement is the cleaner evidence and is what §E.2 cites.

**Every predicted blindness was confirmed by measurement, not asserted.** M8 does
not red reciprocity; M3 does not red the involution row; M1 does not red the
round-trips. Those three are the plan's load-bearing scope statements and they
are now evidence.

---

## A. Existing tests that pin the CURRENT (wrong) periodic behaviour

Inventory from `grep -rn --include='*.py' -E "periodic|Periodic|SpatialWrap|PeriodicWrap" tests/`
(93 hits, 27 files), triaged. Verdict vocabulary per
`.claude/rules/coding-standards.md`: **REWIRE** (behavioural contract that
survives, re-posed), **DELETE** (API-smoke that the successor subsumes),
**KEEP** (characterization / unaffected).

### A.1 The two rows that are RED right now — both REWIRE

| # | `file:line` | What it pins | Verdict |
|---|---|---|---|
| **A-1** | `tests/geometry/test_boundary.py:512-531` `test_periodic_bc_returns_input_unchanged` | "PeriodicBoundary is identity on the angular axis (smoke test) … the spatial pushforward is the caller's responsibility". Realizes through `_realize_for_sn` (`:74`), i.e. `SNMethodSpace.minimal`. | **REWIRE.** The docstring's premise ("the caller's responsibility") is exactly what B3.4c deletes. Re-pose on a face-ful space + the PARTNER's `Γ₊`, and keep the aliasing leg (`psi_in[0,0] = 1e9` must not touch `psi_out`) — that leg is a real contract and is otherwise unpinned once `test_periodic_wrap_operator.py` is trimmed. |
| **A-2** | `tests/sn/operators/test_sn_boundary_realizer.py:695-721` `TestRealizePeriodic::test_periodic_returns_tensor_product_and_passes_through` | The realized TYPE is a 2-factor `TensorProductOperator(PeriodicWrapOperator, IdentityOperator)` AND `op.apply(psi) == psi` on a FULL-FACE `(quad.N, 5)` probe. | **REWIRE, and the value leg must change meaning.** The type leg survives verbatim. The value leg is the pre-B3.4c endomorphism claim — after the narrowing it must be `apply(γ₊ψ|partner) == γ₋ψ|face` on a `(|Γ₊|, 5)` probe, and the sibling negative ("a full-face input is refused / does not emit N rows") belongs here too. |

⚠ `test_boundary.py:63-69` already documents that `_realize_for_sn` **retires
with B3.4c**: *"This helper retires with B3.4c, when periodic narrows and
nothing is left that a faceless space can realize."* A-1 is its last consumer
for periodic; grep `_realize_for_sn` before deleting — it has other callers.

### A.2 Rows that assert the identity BODY of `PeriodicWrapOperator`

| # | `file:line` | What it pins | Verdict |
|---|---|---|---|
| **A-3** | `tests/numerics/test_periodic_wrap_operator.py:27-33` `test_apply_is_identity_various_shapes` | `apply(x) == x` on `(5,)`, `(3,4)`, `(2,3,5)` | **KEEP.** This is the numerics-layer primitive, deliberately shape-agnostic. B3.4c does NOT change the body (§B). |
| **A-4** | `:37-42` `test_apply_transpose_is_identity` | `apply_transpose(x) == x` | **KEEP** — same reason. But see §C6: this row is now the *whole* transpose story for the leaf, and it cannot see the cross-face scatter, which lives in `SNBoundaryOperator`. |
| **A-5** | `:46-50` `test_adjointable_not_invertible_predicates` | `is_adjointable and not is_invertible` | **KEEP.** Unchanged by B3.4c: the leaf is still a self-adjoint identity; what flips is `SpatialWrap.is_adjointable`, a different object. |
| **A-6** | `:54-61` `test_compose_with_identity` | `(P @ I).apply(x) == x` | **KEEP** (algebra closure, orthogonal). |
| **A-7** | `:65-89` `test_apply_returns_fresh_copy` | the safe-aliasing contract | **KEEP**, and it is the reason A-1's aliasing leg is optional rather than mandatory. ⚠ Its docstring cross-references `tests/geometry/test_boundary.py::test_periodic_bc_returns_input_unchanged` (`:75`) — that pointer must be updated when A-1 is re-posed. |
| **A-8** | `tests/numerics/test_operator_capability_predicates.py:192` `("PeriodicWrap", …, False, True, STRUCTURAL_ABSENT)` | not invertible / adjointable / no `inverse()` | **KEEP.** Unchanged. |
| **A-9** | `tests/numerics/test_periodic_wrap_operator.py:1-14` module docstring | claims `apply` "returns the input by reference" — contradicted by A-7 in the same file | **REWIRE (doc).** Pre-existing falsehood, cheap to fix while the file is open; also carries the dangling `PeriodicBoundary.apply` xref (falsehood #3, `:12` and `:68`). |

### A.3 The already-red B3.4c gates (the todo list) — REWIRE per §B

| # | `file:line` | Verdict |
|---|---|---|
| **A-10a** | `tests/geometry/test_bc_equivalence_snapshot.py:773-790` `TestPeriodicLebedev17Snapshot::test_delivers_the_partner_faces_outflow` (strict xfail) | **REWIRE — see §B. The documented edit is a lie; do not apply it verbatim.** |
| **A-10b** | `tests/geometry/test_bc_equivalence_snapshot.py:755-771` `test_the_identity_body_is_correct_on_the_partner_half_trace` | **KEEP** as the body pin. After B3.4c its name should say so (`…_the_realized_wrap_is_the_identity_on_the_partner_half_trace`) — otherwise it and A-10a are two names for one assertion. |
| **A-11** | `tests/sn/operators/test_b3_domain_narrowing.py:312` `_LAWS["periodic"] = (PeriodicBoundary(), True)` + `_B34_XFAIL` `:258-269` | **REWIRE.** The flip is `True → False` + delete `_B34_XFAIL`. This is the campaign's canonical "the marker set IS the todo list" row and it flips honestly: after step 4 the realized law must emit `|Γ₋|` rows for a `|Γ₊|` input *and* not emit `N` rows for a full-face input. ⚠ Leg A of that gate is fed **this face's** `Γ₊` (`np.ones((n_out, 3))`, `:366`) — an all-ones probe, so it is BLIND to which face supplied the data. Leg A stays a shape claim; §C4 carries the value claim. |

### A.4 Rows that MENTION periodic and are unaffected — KEEP, no edit

| # | `file:line` | Why unaffected |
|---|---|---|
| A-12 | `tests/geometry/test_boundary_factors.py:83` `(PeriodicBoundary(), SpatialWrap, ScalarResponse, 1.0, "periodic")` | the (G, R) factor table — B3.4c changes neither factor's TYPE. `:148`'s `pytest.skip` for laws with no `albedo` field still applies. |
| A-13 | `tests/geometry/test_reemission_closure.py:529` `("periodic", PeriodicBoundary(), False)` and `:1138` | "exactly one of G, R is non-trivial" + "not an angular closure". Both survive. |
| A-14 | `tests/geometry/test_bc_universal_invariants.py:709-721` `test_homogeneous_laws_certify_masklessly` | `assert_source_lives_on_incoming_trace` with `q ≡ 0` — law-layer, no realization. |
| A-15 | `tests/geometry/test_boundary_trace_law.py:268`, `test_boundary.py:687,861` | registry-key inventory. |
| A-16 | `tests/geometry/test_bc_errors.py:127-150` | the DIFFUSION refusal (`DiffusionMethodSpace.minimal`), keyed on `isinstance(law.geometry_map, SpatialWrap)` — unchanged and, per explorer §1.4, the free oracle for the whole carve. |
| A-17 | `tests/diffusion/test_boundary_realizer.py:300-305`, `:394` | same refusal + a `LawSum` arm. |
| A-18 | `tests/sn/sweep/test_sweep_acyclicity.py:52` `("periodic","vacuum",False)` | the trace-digraph SCC model; §11.3 defers wiring to its own phase. |
| A-19 | `tests/sn/operators/test_boundary_conditions.py:97-103`, `test_snmesh_realizer_wiring.py:451-468` | assert `BC("periodic")` RAISES + `set(REGISTRY) == {"vacuum","reflective"}`. **These are the gate on finding 0.1(1)** — if B3.4c admits the tag they must be migrated in the SAME commit; if not, they are the proof it did not. |
| A-20 | `tests/sn/operators/test_operator_block_role.py:189` `"periodic": PeriodicBoundary()` | `block_role is BOUNDARY`. The space is face-ful (`:136-146`), so it survives the new guard. ⚠ `_FACE` is `"xmin"`; a `PeriodicBoundary()` defaults to `axis="x"`, so the new axis guard passes — but the file's own note about a latent `xmin`/`xmax` orientation mismatch (carried in agent memory) should be re-checked while it is open. |
| A-21 | `tests/mc/*`, `tests/sn/operators/test_psi_half_coupling.py:756`, `test_fission_operator.py:362`, `test_angular_average_operator.py:19` | prose / MC-side periodic tracking / `NotImplementedError` inventories. No SN realization. |

### A.5 The DECAYED pair — repair, do not flip

| # | `file:line` | Verdict |
|---|---|---|
| **A-22** | `tests/geometry/test_boundary.py:621` `test_wave0_sum_of_realized_bcs_acts_as_weighted_sum` and its sibling `test_operator_sum_of_bcs_acts_as_weighted_sum`, sharing `_MIXED_LAW_XFAIL` (`:605-618`) | **REWIRE — repair, and do NOT count them as B3.4c flips.** |

`[M]` They do not fail for their documented reason:

```
$ .venv/bin/python -O -m pytest tests/geometry/test_boundary.py::test_wave0_sum_of_realized_bcs_acts_as_weighted_sum --runxfail -q
>       white_realized = _realize_for_sn(white, quad)
tests/geometry/test_boundary.py:640:
tests/geometry/test_boundary.py:74: in _realize_for_sn
    return realizer.realize(bc, method_space)
orpheus/sn/boundary/realizer.py:648: in realize
    gamma_out = _outflow_restriction(method_space, "white")
method_space = SNMethodSpace(quadrature=..., face=None, inflow_indices=None, mesh=None, trace=None, outflow_indices=None)
```

The reason text asserts a shape mismatch `(4,2)` vs `(8,2)`; the actual failure
is a faceless method space in the test's own helper. Since white narrowed at
B3.4a the composition IS well-typed — replace `_realize_for_sn` with
`face_method_space(...)` in both bodies, delete `_MIXED_LAW_XFAIL`, and confirm
they go GREEN. If they do, B3.4c has closed two mis-attributed markers that were
falsely charged to it.

---

## B. The flip mechanics of the already-red gate — VERIFIED, and §16.4 is half right

### B.1 `[M]` The gate reds for its documented reason (see also §F)

```
$ .venv/bin/python -O -m pytest tests/geometry/test_bc_equivalence_snapshot.py --runxfail -q
>       np.testing.assert_array_equal(actual, snapshot["psi_in"])
E       AssertionError:
E       Arrays are not equal
E       Mismatched elements: 735 / 735 (100%)
tests/geometry/test_bc_equivalence_snapshot.py:790: AssertionError
1 failed, 43 passed, 1 warning in 1.97s
```

### B.2 §16.4's claim "it will NOT flip by itself" — CONFIRMED BY MEASUREMENT

⭐ **This is no longer a prediction.** Step 4 landed during this session (§0.2)
and the row was re-audited afterwards. `[M]` it still reds, byte-identically:

```
$ .venv/bin/python -O -m pytest tests/geometry/test_bc_equivalence_snapshot.py --runxfail -q
   (AFTER _face_domains + _reflect_trace landed)
E       Mismatched elements: 735 / 735 (100%)
E        [0, 0, 0]: 1.3969355235725678 (ACTUAL), 0.8927739862960189 (DESIRED)
tests/geometry/test_bc_equivalence_snapshot.py:790: AssertionError
1 failed, 43 passed, 1 warning in 0.85s
```

Same 735/735, same first-mismatch values as the pre-step-4 run in §F. The
production change B3.4c exists for is IN, and this "B3.4c gate" did not move
one element. The reasoning:

The xfail body (`:784-790`) calls `_realize(case, space).apply(snapshot["psi_out"])`
— it drives the **realized law object directly**, never `SNBoundaryOperator`. Step 4
changes `_reflect_trace`, i.e. *which array the composition hands the law*. The
law's body is unchanged (`realizer.py` still returns
`PeriodicWrapOperator() & IdentityOperator()`; the working-tree diff adds a guard,
not a body). Feeding it `psi_out` (this face's `Γ₊`, |Γ₊| rows) will still echo
those rows, and `psi_in` is the partner's data — so the row stays RED.

### B.3 §16.4's companion claim "a live companion pins the flip is reachable" — CONFIRMED, and it is ALSO the trap

`test_the_identity_body_is_correct_on_the_partner_half_trace` (`:755-771`) is
green today (it is inside the `43 passed`). Its body is:

```python
space = _verified_space(case, snapshot)                       # :769
actual = _realize(case, space).apply(snapshot["psi_out_partner"])
np.testing.assert_array_equal(actual, snapshot["psi_in"])
```

and `_verified_space` (`:183-185`) builds the SAME space the xfail builds
(`face_method_space(case.build_quadrature(), face=case.face, faces=case.faces)`)
plus two index-set assertions. **The ONLY textual difference between the two
tests is `psi_out_partner` vs `psi_out`.**

### B.4 ⛔ The edit that would be a LIE

> `actual = _realize(case, space).apply(snapshot["psi_out_partner"])` + delete
> the `@pytest.mark.xfail`.

This is exactly what `_B34C_XFAIL` (`:688-691`) instructs. Applying it produces
a duplicate of a row that was **already green before B3.4c**. The marker's
deletion would then be pure ceremony: nothing about the edit is sensitive to
whether step 4 landed, so the "XPASS(strict) failure is the point" mechanism
the campaign relies on silently stops working for this row. That is
vv-principles Mode-8 class 4 (MISATTRIBUTED strict xfail) with the failure mode
inverted — not "xfails for the wrong reason" but "flips for no reason".

### B.5 ✅ The honest edit — re-pose the row at the COMPOSITION, where step 4 lives

The requirement `γ₋ψ|_f = γ₊ψ|_{f'}` is a statement about the **composition**
`ι₋ ∘ law ∘ γ₊`, not about the law. So the row moves to the object that owns
that composition. Concretely, replace the xfail body with:

```python
    def test_the_composition_supplies_the_partner_half_trace(self) -> None:
        r"""**B3.4c** — the REQUIREMENT, at the layer that owns the choice.

        The realized wrap's body is the identity and always was (the sibling
        row above pins it). What B3.4c changes is which half-trace the
        composition ``ι₋ ∘ law ∘ γ_source`` restricts on the way IN. So the
        claim is stated against the law's declared DOMAIN FACE — the one datum
        the composition reads — and then against the value that datum selects.

        Structured so exactly one statement can fail (vv Mode-8 class 4): the
        space cross-check, the pairing premise and the probe activation are the
        LIVE rows above.
        """
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = _verified_space(case, snapshot)
        law = PeriodicBoundary(axis="x")

        # 1. The law's domain is the PARTNER face, derived not asserted.
        assert law.geometry_map.domain_face(case.face) == str(snapshot["partner_face"])

        # 2. Fed the half-trace THAT face names, the realization reproduces
        #    the reference. (Not `x == x`: a scale/permute/average fails.)
        source_face = law.geometry_map.domain_face(case.face)
        probe = {
            case.face: snapshot["psi_out"],
            source_face: snapshot["psi_out_partner"],
        }[source_face]
        actual = SNBoundaryRealizer().realize(law, space).apply(probe)
        np.testing.assert_array_equal(actual, snapshot["psi_in"])
```

Why this is not a lie: leg 1 fails TODAY (`domain_face` did not exist before
this phase; before step 3 the answer would have been `case.face`), and the
`probe` dict is *indexed by the production answer*, so if `domain_face`
regresses to `case.face` the probe becomes `psi_out` and leg 2 reds with the
same 100 % mismatch §F shows. The row is therefore sensitive to the production
change in BOTH directions.

**And the composite-level requirement — the thing step 4 actually builds — is a
NEW test, not an edit of this one.** It cannot live in
`test_bc_equivalence_snapshot.py` at all: that file realizes a law against a
hand-built `SNMethodSpace` and never constructs an `SNBoundaryOperator`. It goes
in `tests/sn/operators/test_sn_boundary_operator.py` — §C4.

### B.6 What would make even the §B.5 edit a lie

If `PeriodicBoundary(axis="x")` in leg 1 is replaced by `case.compose(...)`'s
law or if `source_face` is hard-coded to `"xmax"`. Hard-coding the partner turns
leg 2 back into the already-green sibling and leg 1 into `"xmax" == "xmax"`.
The dict-indexed-by-the-production-answer shape is the load-bearing detail.

---

## C. New gates B3.4c owes

### C0. `[M]` The newly-minted bijection is currently UNGATED

```
$ grep -rn --include='*.py' -E "face_name|face_normal|face_opposite|FACE_NAMES" tests/
… (only `FaceLabel.face_name` / `trace.face_names` hits, plus ONE) …
tests/numerics/test_angular_trace_space.py:469:    from orpheus.numerics.face_layout import face_normal
```

`face_name` (the function), `face_opposite` and `FACE_NAMES` have **zero test
coverage**. The docstring `>>>` examples do NOT run: `pyproject.toml:80` sets
`addopts = "--import-mode=importlib"` with no `--doctest-modules`.

### C1. The cleanup is bit-identical — the EXHAUSTIVE round-trip

New file `tests/numerics/test_face_name_bijection.py`, `pytestmark = [pytest.mark.foundation]`.

```python
import itertools
import pytest
from orpheus.numerics.face_layout import (
    AXIS_NAMES, FACE_NAMES, face_name, face_normal, face_opposite,
)

_AXES = range(len(AXIS_NAMES))
_SIGNS = (-1, +1)


@pytest.mark.parametrize("axis,sign", list(itertools.product(_AXES, _SIGNS)))
def test_render_then_parse_is_the_identity(axis: int, sign: int) -> None:
    """EXHAUSTIVE over the whole bijection — 6 faces, not a sample.

    A convention with a single home is only single-sourced if BOTH directions
    are pinned; a one-way test admits an inverse that disagrees with its own
    forward on one entry.
    """
    assert face_normal(face_name(axis, sign)) == (axis, sign)


@pytest.mark.parametrize("face", FACE_NAMES)
def test_parse_then_render_is_the_identity(face: str) -> None:
    assert face_name(*face_normal(face)) == face


def test_the_inventory_is_the_whole_product() -> None:
    """`FACE_NAMES` is exactly axes x endpoints — no gaps, no duplicates.

    The pre-C5.3 defect this replaces was a hand-listed 4-face table that
    silently lacked z; a set comparison against the derived product is what
    makes that unspellable.
    """
    assert set(FACE_NAMES) == {
        f"{a}{s}" for a in AXIS_NAMES for s in ("min", "max")
    }
    assert len(FACE_NAMES) == len(set(FACE_NAMES)) == 2 * len(AXIS_NAMES)
```

`[M]` all six rows currently hold:

```
FACE_NAMES = ('xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax')
name->normal->name round-trip failures: []
normal->name round-trip failures: []
```

**Plus the bit-identity leg for the repoint** — the reason the cleanup is
claimed no-op. Three of the five sites rendered
`f"{axis}{'max' if outward_sign == +1 else 'min'}"`; pin the successor against
that retired expression *transcribed in the test*, exhaustively:

```python
@pytest.mark.parametrize("axis,sign", list(itertools.product(_AXES, _SIGNS)))
def test_the_bijection_reproduces_the_five_retired_transcriptions(axis, sign):
    """The repoint is BIT-IDENTICAL — the retired expression, transcribed.

    A refactor that collapses N transcriptions of a convention onto one
    primitive is only a no-op if the primitive agrees with all N. The retired
    render is four tokens; restating it here is cheaper than trusting that no
    site disagreed, and it is a DIFFERENT expression from the successor's
    dict lookup (structural, not procedural, independence).
    """
    legacy = f"{AXIS_NAMES[axis]}{'max' if sign == +1 else 'min'}"
    assert face_name(axis, sign) == legacy
    assert face_normal(legacy) == (axis, sign)          # the .endswith("max") twin
```

**Negative legs** (`[M]` all already raise `ValueError`):

```
face_name(3, 1)  -> ValueError: axis index 3 is outside the named-axis inventory
face_name(0, 0)  -> ValueError: outward_sign must be +1 ... or -1
face_name(0, 2)  -> ValueError
face_normal('wmin'/'xmid'/'x'/''/'minx'/'xminmax') -> ValueError  (all six)
```

`'minx'` and `'xminmax'` matter: `face_normal` is implemented with `endswith`,
so a suffix-shaped name is the plausible-wrong acceptance.

### C2. `face_opposite` is an involution and axis-preserving

```python
@pytest.mark.parametrize("face", FACE_NAMES)
def test_opposite_is_an_involution(face: str) -> None:
    assert face_opposite(face_opposite(face)) == face


@pytest.mark.parametrize("face", FACE_NAMES)
def test_opposite_preserves_the_axis_and_flips_the_normal(face: str) -> None:
    """Both halves, because either alone admits a wrong map.

    Axis-only would admit the identity; sign-only would admit xmin -> ymax.
    `[M]` |G-(xmin)| == |G+(ymax)| == 12 on level_symmetric(4), so a
    cross-axis partner is SHAPE-LEGAL and a shape assertion cannot see it
    (vv Mode 12) — this is where it is caught, by NAME, before any array.
    """
    axis, sign = face_normal(face)
    assert face_normal(face_opposite(face)) == (axis, -sign)


@pytest.mark.parametrize("face", FACE_NAMES)
def test_opposite_has_no_fixed_point(face: str) -> None:
    """A free action — the torus, not the orbifold.

    This is the one line that separates `SpatialWrap` from `SpecularMirror`
    at the level of the map: the mirror's quotient HAS fixed points (the
    plane), the wrap's does not.
    """
    assert face_opposite(face) != face
```

`[M]` `involution failures: []`, `axis-preservation failures: []`,
`sign-flip failures: []`, `no fixed points: []`.

### C2b. ⚠ The SIXTH transcription the cleanup did not collapse

`orpheus/transport/mesh/axis.py:147` still renders
`f"{AXIS_NAMES[self.axis_index]}{suffix}"` from a THREE-valued endpoint
(`min`/`max`/`outer`). It is not a duplicate of `face_name` (its domain is
larger), but the *render* is duplicated. `[M]` they agree today:

```
FaceLabel(0,'min').face_name = 'xmin'   face_name = 'xmin'   agree=True
FaceLabel(0,'max').face_name = 'xmax'   face_name = 'xmax'   agree=True
FaceLabel(1,'min') … (2,'max')          all agree=True
FaceLabel(0,'outer').face_name = xmax   <- three-valued endpoint, no sign counterpart
```

Whether `FaceLabel.face_name` should delegate is a design call for the main
agent. **Either way** the agreement wants a gate, because it is the seam every
mesh-built face name crosses:

```python
@pytest.mark.parametrize("axis,ep,sign", [(a, ep, s) for a in _AXES
                                          for ep, s in (("min", -1), ("max", +1))])
def test_facelabel_render_agrees_with_the_bijection(axis, ep, sign):
    """The mesh-side crosswalk and the numerics-side bijection are ONE world.

    `FaceLabel` carries a THREE-valued endpoint ("outer" renders as max), so
    it is not simply `face_name`; on the two-valued subset they must agree or
    a mesh-built face name and a law-derived partner name refer to different
    faces. Pinned on the subset where the question is even askable.
    """
    from orpheus.transport.mesh.axis import FaceLabel
    assert FaceLabel(axis, ep).face_name == face_name(axis, sign)
```

### C3. The `domain_face` member — default identity, `SpatialWrap` RAISES on a mismatched axis

New file `tests/geometry/test_boundary_domain_face.py`, `foundation`.

```python
_MAPS_THAT_STAY = [
    pytest.param(IdentityMap(), id="identity"),
    pytest.param(SpecularMirror(axis="x"), id="specular-x"),
    pytest.param(SpecularMirror(axis="y"), id="specular-y"),
]


@pytest.mark.parametrize("gmap", _MAPS_THAT_STAY)
@pytest.mark.parametrize("face", FACE_NAMES)
def test_a_non_quotient_map_consumes_the_face_it_is_installed_on(gmap, face):
    r"""Every map but the wrap answers `face` — and it must for EVERY face.

    Parametrised over all six faces including the mismatched ones
    (`SpecularMirror(axis="y")` on `xmin`): the mirror deliberately does NOT
    reconcile its axis here — that check lives at realization, where the
    reflection TABLE can diagnose it — so a name-only refusal added here
    would be a second, weaker authority that can disagree with the first.
    """
    assert gmap.domain_face(face) == face


@pytest.mark.parametrize("axis", AXIS_NAMES)
def test_the_wrap_consumes_the_OPPOSITE_face(axis):
    for sign in (-1, +1):
        face = face_name(AXIS_NAMES.index(axis), sign)
        assert SpatialWrap(axis=axis).domain_face(face) == face_opposite(face)
        # ...and it is genuinely a DIFFERENT face — the one property that
        # separates a quotient from a wall.
        assert SpatialWrap(axis=axis).domain_face(face) != face


@pytest.mark.parametrize("wrap_axis,face", [
    ("y", "xmin"), ("y", "xmax"), ("x", "ymin"), ("x", "ymax"),
    ("z", "xmin"), ("x", "zmax"),
])
def test_a_wrap_on_a_foreign_axis_REFUSES(wrap_axis, face):
    """ERR-041 class: two encodings of one orientation that disagree.

    A wrap along y installed on an x face identifies nothing. Returning
    `face` would realize it as a bare identity on the wrong half-trace and
    report nothing — the silent shape B3.4a fixed for white.
    """
    with pytest.raises(BoundaryError) as exc:
        SpatialWrap(axis=wrap_axis).domain_face(face)
    assert exc.value.law == "periodic"
    assert wrap_axis in str(exc.value) and face in str(exc.value)


def test_an_unparseable_face_REFUSES_as_a_BoundaryError():
    """Not a bare ValueError: the caller is a realizer, and the campaign's
    contract is that a boundary refusal names its law."""
    with pytest.raises(BoundaryError):
        SpatialWrap(axis="x").domain_face("bogus")
```

`[M]` the axis guard already fires with the right message:

```
realize(PeriodicBoundary(axis='y'), 'xmin')
 -> BoundaryError: SpatialWrap declares axis='y' but is installed on face 'xmin', which lies
    on axis 'x'. A wrap identifies the two faces NORMAL to its own axis; ...
```

**Protocol-completeness leg** (the illegal-states gate for step 3):

```python
def test_every_geometry_map_answers_domain_face():
    """A map that cannot name its own domain face is not realizable.

    The Protocol gained a THIRD member; a concrete that forgot it fails at
    the first realization, not at definition. This is the row that makes a
    future map's omission loud.
    """
    for gmap in (IdentityMap(), SpecularMirror(), SpatialWrap()):
        assert isinstance(gmap.domain_face("xmax"), str)
```

### C4. ⭐ The end-to-end physics — `B·ψ` crosses the faces

**This is the gate for step 4, the one thing not yet built, and it does not
exist anywhere today.** Home: `tests/sn/operators/test_sn_boundary_operator.py`,
next to the RG-1..RG-4 family (`TestNarrowedLawDomain`, `:479`).

`[M]` the fixture is buildable today and RED today:

```
faces: ('xmin', 'xmax')
  xmin: Gamma- = [4, 5, 6, 7]  Gamma+ = [0, 1, 2, 3]
  xmax: Gamma- = [0, 1, 2, 3]  Gamma+ = [4, 5, 6, 7]
injected kinds: {'xmin': 'periodic', 'xmax': 'periodic'}

face xmin: emitted rows on Gamma-(xmin) = [4, 5, 6, 7]
  == PARTNER's Gamma+ data? False
  == OWN     Gamma+ data? True
  max|got - partner| = 1.566388e+00
  max|got - own    | = 0.000000e+00
  off-Gamma- rows all zero? True
face xmax:  (identical verdicts)
```

```python
def _periodic_slab(n_ordinates: int = 8) -> SNMesh:
    r"""A slab with periodic INJECTED on both faces.

    ⚠ `BC("periodic")` cannot build this: `SNMesh.BOUNDARY_OPERATOR_REGISTRY`
    is `{"vacuum", "reflective"}` (augmented_mesh.py:175-178) and admitting
    periodic is #189. So the law is installed through the same
    `realize_boundary_law` hook the tag path uses — one layer below the
    registry, and the ONLY layer at which this claim is stateable until #189.
    `gauss_legendre` because it is a DISCRIMINATING quadrature (§D).
    """
    quad = Quadrature.gauss_legendre(n_ordinates=n_ordinates)
    mesh = Mesh1D(edges=np.linspace(0.0, 2.0, 5), mat_ids=np.zeros(4, dtype=int))
    sn = SNMesh(mesh, quad, placeholder_materials())
    for face in ("xmin", "xmax"):
        sn.bc[face] = sn.realize_boundary_law(PeriodicBoundary(axis="x"), face)
    return sn


def _independently_seeded_trace(sn: SNMesh, seed: int) -> TimedFullField:
    r"""A boundary state whose two faces carry DIFFERENT data.

    §16.4's central lesson, one level up: with one shared draw, "the partner's
    outflow" and "this face's outflow" are different rows of one array and a
    per-face endomorphism looks defensible. The independence is what makes the
    claim observable — and it is ASSERTED below, not assumed.
    """
    z = TimedFullField.zeros(interior=AngularFlux,
                             boundary=AngularBoundaryFlux, mesh=sn)
    rng = np.random.default_rng(seed)
    return replace(z, boundary=replace(
        z.boundary, values=rng.uniform(0.5, 2.0, size=z.boundary.values.shape)))


class TestPeriodicIsACrossFaceCoupling:
    r"""**B3.4c** — ``B`` is block-STRUCTURED, not block-diagonal, over faces."""

    def test_the_two_faces_carry_different_data(self) -> None:
        """ACTIVATION, outside the claim (vv Mode-8 class 4 discipline).

        If the two faces ever coincided, every row below would be green for a
        per-face endomorphism and the class would assert nothing.
        """
        sn = _periodic_slab()
        psi = _independently_seeded_trace(sn, seed=20260801)
        tr = sn.angular_trace
        assert not np.array_equal(
            psi.boundary.face_view("xmin")[tr.outflow_indices_for_face("xmin")],
            psi.boundary.face_view("xmax")[tr.outflow_indices_for_face("xmax")],
        )

    @pytest.mark.parametrize("face,partner", [("xmin", "xmax"), ("xmax", "xmin")])
    def test_the_inflow_is_the_PARTNERS_outflow(self, face, partner) -> None:
        r"""``(Bψ)|_{Γ₋(f)} == ψ|_{Γ₊(f')}`` — the whole of B3.4c.

        BOTH directions are parametrised: a one-face gate would pass on an
        implementation that wrapped xmin from xmax and left xmax echoing
        itself, which is exactly what a `face_opposite` used in one branch
        and forgotten in the other produces.

        `assert_array_equal`: a translation is a gather. Any tolerance here
        would admit the failure the row exists to catch.
        """
        sn = _periodic_slab()
        psi = _independently_seeded_trace(sn, seed=20260801)
        tr = sn.angular_trace
        got = SNBoundaryOperator(sn).apply(psi).boundary.face_view(face)[
            tr.inflow_indices_for_face(face)]
        want = psi.boundary.face_view(partner)[
            tr.outflow_indices_for_face(partner)]
        np.testing.assert_array_equal(got, want)

    @pytest.mark.parametrize("face", ["xmin", "xmax"])
    def test_it_is_NOT_the_faces_own_outflow(self, face) -> None:
        r"""The anti-Mode-12 leg, and the one that reds TODAY.

        `[M]` `Γ₋(xmin) == Γ₊(xmax)` as index SETS on every quadrature in the
        tree, so the emitted ROWS are the same either way and only the VALUES
        discriminate. `max|got − own| = 0.0` today, `max|got − partner| = 1.57`.
        """
        sn = _periodic_slab()
        psi = _independently_seeded_trace(sn, seed=20260801)
        tr = sn.angular_trace
        got = SNBoundaryOperator(sn).apply(psi).boundary.face_view(face)[
            tr.inflow_indices_for_face(face)]
        own = psi.boundary.face_view(face)[tr.outflow_indices_for_face(face)]
        assert not np.array_equal(got, own), (
            f"{face}: B echoed this face's own Γ₊ — the pre-B3.4c "
            f"per-face endomorphism. |Γ₊| == |Γ₋| here, so no shape check "
            f"can see this (vv Mode 12)."
        )

    @pytest.mark.parametrize("face", ["xmin", "xmax"])
    def test_nothing_lands_off_gamma_minus(self, face) -> None:
        """The CODOMAIN half of the narrowing — a wrap that scattered onto
        Γ₊ would corrupt the outflow-definition residual, which carries no B
        term at all (`boundary.py:326-352`)."""
        sn = _periodic_slab()
        psi = _independently_seeded_trace(sn, seed=20260801)
        tr = sn.angular_trace
        out = SNBoundaryOperator(sn).apply(psi).boundary.face_view(face)
        off = np.setdiff1d(np.arange(out.shape[0]),
                           tr.inflow_indices_for_face(face))
        np.testing.assert_array_equal(out[off], np.zeros_like(out[off]))
```

**2-D companion — MANDATORY, and it is the only place the axis choice is
observable.** `[M]` on `level_symmetric(4)` a wrong-AXIS partner is shape-legal:

```
faces: ('xmin', 'xmax', 'ymin', 'ymax')
kinds: {'xmin': 'periodic', 'xmax': 'periodic', 'ymin': 'reflective', 'ymax': 'reflective'}
xmin: shapes got=(12, 1, 3) partner=(12, 1, 3) self=(12, 1, 3)
   == partner? False   == self? True     max|got-partner| = 1.8423e+00
|G-(xmin)| = 12  |G+(ymax)| = 12  -> a wrong-axis partner would be SHAPE-LEGAL: True
G-(xmin) == G+(ymax) as SETS? False
```

So the 2-D row must (a) install periodic on the x-pair only, (b) assert the
y-pair's reflective images are BYTE-IDENTICAL to a control run with vacuum on
x (the "periodic did not leak across axes" leg), and (c) assert
`Γ₋(xmin) != Γ₊(ymax)` as index sets so the reader sees why a shape check is
insufficient.

### C5. The transpose — `⟨Bx, y⟩ == ⟨x, Bᵀy⟩`, x and y independently seeded

`[M]` today's operator is ALREADY perfectly reciprocal, because forward and
transpose are wrong in the same way:

```
  <Bx,y> = 7.706205854030817
  <x,B'y>= 7.706205854030817
  |diff| = 0.000000e+00   rel = 0.000e+00
```

**So reciprocity ALONE is blind to the B3.4c defect** — it is a
forward/transpose *consistency* check, not a correctness check. It becomes a
real catcher only for the transpose-side mutations. `[M]` synthetic
whole-trace matrices, x and y independently seeded:

```
gauss_legendre(8)  N=8  |G-|=4
  TODAY   : F_today   with T_same_face   |gap|=8.882e-16  rel=1.153e-16
  B3.4c   : F_correct with T_correct     |gap|=8.882e-16  rel=1.205e-16
  MUTANT-1: F_correct with T_same_face   |gap|=3.331e-01  rel=4.322e-02
  MUTANT-2: F_correct with T_wrong_half  |gap|=1.725e+00  rel=2.340e-01
  MUTANT-3: F_today   with T_correct     |gap|=3.331e-01  rel=4.322e-02
lebedev(17)  N=110  |G-|=49
  TODAY   : F_today   with T_same_face   |gap|=1.421e-14  rel=1.454e-16
  B3.4c   : F_correct with T_correct     |gap|=0.000e+00  rel=0.000e+00
  MUTANT-1: F_correct with T_same_face   |gap|=1.492e+00  rel=1.526e-02
  MUTANT-2: F_correct with T_wrong_half  |gap|=1.126e+01  rel=1.048e-01
  MUTANT-3: F_today   with T_correct     |gap|=1.492e+00  rel=1.526e-02
```

```python
@pytest.mark.parametrize("n_ordinates", [4, 8])
def test_euclidean_reciprocity_on_a_periodic_slab(n_ordinates) -> None:
    r"""``<Bx, y> == <x, Bᵀy>`` with x and y INDEPENDENTLY seeded.

    Scope, stated because it is not obvious: this row is BLIND to the B3.4c
    defect itself — `[M]` the pre-B3.4c pairing (own-face forward, own-face
    transpose) reciprocates EXACTLY (rel 1.2e-16), because both legs are
    wrong the same way. What it catches is the forward/transpose MISMATCH the
    carve can introduce: `[M]` a transpose that scatters back over this face's
    Γ₊ instead of the partner's reads rel 4.3e-2 (GL8) / 1.5e-2 (lebedev17),
    and one that scatters over the partner's Γ₋ reads 2.3e-1 / 1.0e-1.

    Independent seeds are load-bearing: with x == y the bilinear form is
    symmetric by construction and every mutant above collapses.

    `rtol=0` / `atol = n_terms · eps`: both sides are the SAME sum of the
    SAME products in a different association order; there is no cancellation
    (all draws positive).
    """
    sn = _periodic_slab(n_ordinates=n_ordinates)
    B = SNBoundaryOperator(sn)
    x = _independently_seeded_trace(sn, seed=1)
    y = _independently_seeded_trace(sn, seed=2)
    lhs = float(np.sum(B.apply(x).boundary.values * y.boundary.values))
    rhs = float(np.sum(x.boundary.values * B.apply_transpose(y).boundary.values))
    n_terms = x.boundary.values.size
    np.testing.assert_allclose(lhs, rhs, rtol=0.0,
                               atol=n_terms * float(np.finfo(np.float64).eps)
                                    * max(abs(lhs), abs(rhs)))


def test_the_transpose_scatters_over_the_PARTNERS_outflow() -> None:
    r"""The OBJECT gate the bilinear form cannot give (vv Mode 12).

    Reciprocity is a scalar functional; the whole cross-face structure is
    inside its invariance group whenever forward and transpose agree. So pin
    the transpose's SUPPORT directly: fed a state supported only on Γ₋(xmin),
    Bᵀ must deposit on Γ₊(xmax) — the PARTNER's slot — and leave xmin's own
    Γ₊ at zero.
    """
    sn = _periodic_slab()
    tr = sn.angular_trace
    z = TimedFullField.zeros(interior=AngularFlux,
                             boundary=AngularBoundaryFlux, mesh=sn)
    y = replace(z, boundary=replace(z.boundary, values=z.boundary.values.copy()))
    y.boundary.face_view("xmin")[tr.inflow_indices_for_face("xmin")] = 1.0

    out = SNBoundaryOperator(sn).apply_transpose(y).boundary
    np.testing.assert_array_equal(
        out.face_view("xmax")[tr.outflow_indices_for_face("xmax")],
        np.ones_like(out.face_view("xmax")[tr.outflow_indices_for_face("xmax")]),
    )
    np.testing.assert_array_equal(
        out.face_view("xmin"), np.zeros_like(out.face_view("xmin")),
    )
```

### C6. `SpatialWrap.is_adjointable` `False → True` — where the gate must NOT be

`[M]` the flip is unobservable at the composite:

```
realize_boundary_law(PeriodicBoundary(), 'xmin') -> OK: TensorProductOperator
  law: PeriodicBoundary(axis='x')  kind: periodic
  is_adjointable: True        <-- BEFORE the flip
B.is_adjointable: True        <-- SNBoundaryOperator, periodic on both faces
```

`SNBoundaryOperator.is_adjointable` (`boundary.py:291-310`) intersects each
**realized** law's predicate — never the law FACTOR. So:

```python
def test_the_wrap_declares_an_honest_transpose() -> None:
    r"""**B3.4c** — the factor's `is_adjointable` flipped, and it is a THEOREM.

    Scope, MEASURED: this row is a claim about the LAW FACTOR only. It cannot
    be stated against the realized operator or against `B.is_adjointable`,
    because both read the realized `PeriodicWrapOperator` (identity body) and
    reported `True` even while `SpatialWrap` reported `False` — the two
    sources of truth doc falsehood #4 named. The teeth for the transpose
    itself are `test_the_transpose_scatters_over_the_PARTNERS_outflow`.
    """
    assert SpatialWrap(axis="x").is_adjointable


def test_every_deck_transformation_is_adjointable() -> None:
    """It is a theorem for a measure-preserving bijection, so the tier is
    uniform — and after B3.4c there is no map left reporting an
    implementation gap. A future `False` here is a claim that some deck
    element is NOT a bijection, which needs its own argument."""
    for gmap in (IdentityMap(), SpecularMirror(), SpatialWrap()):
        assert gmap.is_adjointable
```

### C7. What pins that `B` is block-STRUCTURED — and the gates whose REASON decays

**YES, there are existing block-diagonality gates — and the dangerous part is
that they will all stay GREEN.** `[M]` all three are posed on reflective-only
fixtures (`_sn("SLB", (BC.reflective, BC.reflective))`), so B3.4c cannot red
them; what decays is their stated justification. This is vv Mode-8 class 7 (the
DECAYED marker with a half-life): nothing in the test, the tag, or CI notices.

| `file:line` | claim | after B3.4c |
|---|---|---|
| `tests/sn/operators/test_sn_boundary_operator.py:332-355` `test_block_diagonal_no_face_mixing` | "a perturbation on ONE face's input slot affects ONLY that face's output (``B`` is block-diagonal over faces — it never mixes faces)" | **GREEN** (reflective fixture). The claim must be re-scoped: block-diagonality is a property of the installed LAWS, not of ``B``. Rename to `test_a_specular_pair_does_not_mix_faces` or add the law scope to the docstring. |
| `:16` module docstring — *"**block-diagonal over faces** — a single-face perturbation stays on that face"* | a stated property of ``B`` | **FALSE as stated.** Repair with the file. |
| `:656-700` `TestFaceRestrictedReflect` (`test_subset_reflects_only_selected_faces`, `test_subset_partitions_the_whole_trace_reflect`) | *"``B`` is block-diagonal over faces, so the subset MUST be the EXACT restriction … no coupling dropped, no term double-counted"* | **GREEN**, and its REASON is now false in general. ⚠ Load-bearing: this is the Gauss-Seidel octant schedule's exactness argument. |

And in production, the same claim three more times:

* `orpheus/sn/operators/boundary.py:44-45` — *"``B`` is therefore block-diagonal
  over faces — it never mixes faces."*
* `:388-392` — *"``B`` is block-diagonal over faces, so the subset action is the
  EXACT restriction (no cross-face coupling is dropped)."*
* `:619-622` (`split`) — periodic listed among "never-reflected faces" lagged
  into `B_upper`, with the partition argued as exact "by construction".

B3.4c owes a POSITIVE gate that makes the new structure a checked claim:

```python
def test_B_is_block_structured_not_block_diagonal() -> None:
    r"""The structural claim the docstring makes, as an assertion.

    Excite ONLY xmax's Γ₊ and observe xmin's Γ₋ respond. Pre-B3.4c this is
    impossible by construction ("it never mixes faces"); post-B3.4c it is the
    definition of a periodic pair.

    ⚠ ZERO BASE, not a difference of two random states. `[M]` the differencing
    form `B(base + e) − B(base)` is NOT bit-exact — 2/4 entries land at
    0.9999999999999998 / 1.0000000000000002, the catastrophic cancellation of
    `(x+1) − x`. Exploiting linearity instead makes the assertion exact and
    keeps `assert_array_equal`, which is what a gather deserves.

    The third leg is free and pins the DIRECTION: xmax's own Γ₋ reads xmin's
    Γ₊, which is zero here, so the whole xmax slot must stay zero. A wrap
    wired symmetrically-but-wrongly (both faces reading xmax) would light it up.
    """
    sn = _periodic_slab()
    tr = sn.angular_trace
    z = TimedFullField.zeros(interior=AngularFlux,
                             boundary=AngularBoundaryFlux, mesh=sn)
    bumped = replace(z, boundary=replace(z.boundary, values=z.boundary.values.copy()))
    bumped.boundary.face_view("xmax")[tr.outflow_indices_for_face("xmax")] = 1.0

    out = SNBoundaryOperator(sn).apply(bumped).boundary
    sel = out.face_view("xmin")[tr.inflow_indices_for_face("xmin")]
    np.testing.assert_array_equal(sel, np.ones_like(sel))          # [M] bit-exact
    off = np.setdiff1d(np.arange(out.face_view("xmin").shape[0]),
                       tr.inflow_indices_for_face("xmin"))
    np.testing.assert_array_equal(out.face_view("xmin")[off],
                                  np.zeros_like(out.face_view("xmin")[off]))
    np.testing.assert_array_equal(out.face_view("xmax"),
                                  np.zeros_like(out.face_view("xmax")))


def test_a_face_subset_reflect_on_a_WRAP_face(self) -> None:
    r"""HONEST SCOPE — `TestFaceRestrictedReflect`'s reason is now false in
    general, and this row is what stops it being silently re-generalised.

    `[M]` `TestFaceRestrictedReflect` (`test_sn_boundary_operator.py:656-700`)
    is posed on `(BC.reflective, BC.reflective)` and stays GREEN — its
    JUSTIFICATION ("``B`` is block-diagonal over faces, so the subset MUST be
    the EXACT restriction") is what decays. `_reflect_trace(faces={"xmin"})`
    was exact because B never mixed faces; with a wrap installed the xmin block
    READS xmax, so restricting the FACE SET no longer restricts the coupling.

    Written as a documented expectation, NOT a bug. B3.4c must choose:
      (a) the subset action stays exact because the wrap reads the partner's
          INPUT slot (which the subset does not mask) — then assert
          `only_xmin.face_view("xmin") == full.face_view("xmin")` here too and
          re-scope the sibling's docstring from "``B``" to "these laws"; or
      (b) the subset action RAISES for a wrap face — then assert that.
    Either way it is pinned; what must not happen is the sibling's exactness
    argument being read as covering the wrap.
    """
```

⚠ **Scope flag for the plan of record, not for me to resolve:** §11.3 already
says B3.4c "lags periodic into `B_upper`, which is what today's split already
does and is correct-but-not-minimal". That is only true while `split`'s
`upper_rows` computation is face-local. Verify it after step 4 and pin whichever
answer is chosen — an un-gated exactness claim in a splitting is the ERR-056
shape.

**Retirement/rescoping list for C7** (a green gate whose reason decayed is the
Mode-8 class-7 shape — the repair is doc-side, but it is a REQUIRED deliverable):

| site | action |
|---|---|
| `tests/sn/operators/test_sn_boundary_operator.py:16` | re-scope "``B`` is block-diagonal over faces" → "block-diagonal for the laws in `_CASES`; a periodic pair is block-STRUCTURED (see `TestPeriodicIsACrossFaceCoupling`)" |
| `:332-334` `test_block_diagonal_no_face_mixing` | rename + re-scope the docstring to the specular fixture; do NOT delete (the claim is real for those laws) |
| `:662-666` `TestFaceRestrictedReflect` | re-scope the class docstring; add the wrap row above |
| `orpheus/sn/operators/boundary.py:44-45`, `:388-392`, `:619-622` | the production half of the same three claims |

### C8. The missing-partner negative — `[M]` currently SILENT

```
sphere faces: ('xmax',)
face_opposite('xmax') = xmin  -> present on the mesh? False
realize(periodic, 'xmax') on a SPHERE -> ACCEPTED: TensorProductOperator  <-- today: silent
```

A curvilinear mesh has ONE face. A periodic law installed there names a partner
that does not exist, and today's realization accepts it and returns a wrap
whose domain is a face the trace has no slot for. Step 4 will then index a
missing face at apply time — a `KeyError` deep in `_reflect_trace`, or worse,
silence.

```python
def test_a_wrap_whose_partner_face_is_absent_is_REFUSED() -> None:
    r"""A one-faced mesh is not a torus.

    `[M]` today this realizes SILENTLY: a sphere's trace carries only `xmax`,
    `face_opposite("xmax") == "xmin"` is absent, and the realizer returns a
    wrap regardless. Post-B3.4c the composition would look that face up. The
    refusal belongs at REALIZATION, next to the identification guard, so the
    diagnosis names the topology rather than surfacing as a KeyError.
    """
    sn = _sphere()  # Mesh1D(coord=SPHERICAL), faces == ("xmax",)
    with pytest.raises(BoundaryError) as exc:
        sn.realize_boundary_law(PeriodicBoundary(axis="x"), "xmax")
    assert exc.value.law == "periodic"
    assert "xmin" in str(exc.value)
```

### C9. The tangential-band fixture — `product(2,4)`

`[M]` the one quadrature whose tangential cosines are round-off, not exact zero
(§16.6's discriminating fixture):

```
xmin: |G-|=2 [0, 4]  |G+|=2 [2, 6]  tangential=4 [1, 3, 5, 7]
    raw mu_x on tangential = [ 4.99959962e-17 -1.49987989e-16  4.99959962e-17 -1.49987989e-16]
xmax: |G-|=2 [2, 6]  |G+|=2 [0, 4]  tangential=4 [1, 3, 5, 7]
G-(xmin) == G+(xmax)? True
```

So `_assert_wrap_identification` must be exercised on `product(2,4)`, where 4 of
8 ordinates are in NEITHER half-trace. `[M]` the identification survives both
the band and the strict `> 0.0` twin:

```
product(2,4)         band: G-(xmin)==G+(xmax) -> True   strict: True   |G-| band=2 strict=4
gauss_legendre(8)    band: G-(xmin)==G+(xmax) -> True   strict: True   |G-| band=4 strict=4
lebedev(17)          band: G-(xmin)==G+(xmax) -> True   strict: True   |G-| band=49 strict=49
```

**⚠ HONEST SCOPE — a named blind mutation.** The `TANGENTIAL_EPS → 0.0` twin
does NOT red the identification guard on any quadrature: the strict cut changes
BOTH half-traces symmetrically, so the equality is preserved (band `|Γ₋| = 2` →
strict `4` on `product(2,4)`, and still equal). The guard is therefore blind to
the B3.4a classification twin **by construction**, and correctly so — it is a
different claim. What DOES catch it for periodic is
`test_bc_equivalence_snapshot.py::_verified_space`'s index-set cross-check
(`:186-200`) on a quadrature where the band matters, i.e. `product(2,4)` —
which that file does not carry (§16.6). Adding a `periodic_product24` snapshot
case is the only way to close it there; **I recommend NOT adding it** (it would
duplicate §C9's realization gate) and instead stating the scope in the guard's
own docstring, per §16.6's precedent.

---

## D. Mode-12 / config-blindness audit

Per gate, what it is blind to.

| gate | BLIND to | why | covered by |
|---|---|---|---|
| C1 round-trip | a wrong CONVENTION shared by both directions (if `min` meant `+1` everywhere) | a bijection's round-trip is invariant under relabelling both maps | C1's `legacy` transcription leg + `test_face_normals_z_derived_from_axis_names` (`test_angular_trace_space.py:455`), which pins the ABSOLUTE table |
| C2 involution | a partner on the WRONG AXIS if only involution is asserted (`xmin↔ymax` is also an involution) | involution is a property of the permutation, not of the axis | C2's axis-preservation leg; C4's 2-D row |
| C2 no-fixed-point | nothing structural; it is the `SpecularMirror` discriminator | — | — |
| C3 `domain_face` | whether the answer is USED — a `domain_face` that is correct and never read is green | it is a pure name-level claim | C4 (values) + E-M4 |
| C4 value rows | **the tangential band**: on `product(2,4)` a strict-`>0.0` twin changes both halves symmetrically, so the identification survives | the error is inside the equality functional's invariance group | stated as honest scope (§C9); NOT closeable by a value gate |
| C4 value rows | a bug that swaps BOTH faces consistently (xmin reads xmax's Γ₋ instead of Γ₊, and vice versa) — but `[M]` `Γ₋(partner) = Γ₊(face)`, so this is the "own face" mutation already covered | — | `test_it_is_NOT_the_faces_own_outflow` |
| C4 shape legs | which face supplied the data — `[M]` `|Γ₊| == |Γ₋|` on **every** quadrature × face in the tree, and `Γ₋(xmin)` and `Γ₊(xmax)` are the SAME index set, so shape AND emitted-row-index are both invariant | the discriminator has to be the VALUE, on independently-seeded faces | the whole design of C4 |
| C5 reciprocity | **the B3.4c defect itself** — `[M]` rel 1.15e-16 on the pre-B3.4c pairing | forward and transpose are wrong in the same way; the bilinear form's invariance group contains the shared error | C4 (forward) + `test_the_transpose_scatters_over_the_PARTNERS_outflow` (support, an OBJECT gate) |
| C5 reciprocity | a transpose wrong by a *symmetric* relabelling within the correct face | any `T` with `T = Fᵀ` reciprocates by definition | the support gate |
| C6 factor flip | everything downstream — `[M]` the realized operator and `B.is_adjointable` already answered `True` | the composite predicate reads the realized operator, never the factor | declared as scope in the test's own docstring |
| C7 block-structure | a leak in the OTHER direction (xmax → xmin only) | the perturbation is one-sided | parametrise both directions, as C4 does |

### D.1 `[M]` Which quadratures DISCRIMINATE

`§16.5`'s hazard, re-measured. `mirror(Γ₋)` vs `sorted(Γ₊)`:

```
quadrature           face   identity == mirror-on-Gamma-?
gauss_legendre(8)    xmin   False    mirror(G-)=[3, 2, 1, 0] sorted(G+)=[0, 1, 2, 3]
gauss_legendre(8)    xmax   False    mirror(G-)=[7, 6, 5, 4] sorted(G+)=[4, 5, 6, 7]
gauss_legendre(4)    xmin   False    mirror(G-)=[1, 0] sorted(G+)=[0, 1]
gauss_legendre(4)    xmax   False    mirror(G-)=[3, 2] sorted(G+)=[2, 3]
product(2,4)         xmin   True     mirror(G-)=[2, 6] sorted(G+)=[2, 6]
product(2,4)         xmax   True     mirror(G-)=[0, 4] sorted(G+)=[0, 4]
level_symmetric(4)   xmin   True     (12 entries, identical)
level_symmetric(4)   xmax   True     (12 entries, identical)
level_symmetric(6)   xmin   True     (24 entries, identical)
level_symmetric(6)   xmax   True     (24 entries, identical)
lebedev(17)          xmin   False    (49 entries; positions 9/21/33/36 permuted)
lebedev(17)          xmax   False    (49 entries; positions 9/21/33/36 permuted)
```

**DISCRIMINATING: `gauss_legendre(4)`, `gauss_legendre(8)` (the mirror REVERSES
order), `lebedev(17)` (the mirror SCRAMBLES).
BLIND: `product(2,4)`, `level_symmetric(4)`, `level_symmetric(6)`.**

Every value gate above uses `gauss_legendre` or `lebedev`. `product(2,4)`
appears exactly once (§C9) and only for the tangential band, where it is the
discriminating fixture.

### D.2 A blindness that does NOT apply, and why it is worth recording

The obvious §16.5 worry for periodic is "did someone apply the mirror as well
as the wrap?" `[M]` that mutation is **structurally unrealizable**: the mirror
of `Γ₋(xmin) = {4,5,6,7}` is `{3,2,1,0}`, which is not inside
`Γ₊(xmax) = {4,5,6,7}`. The probe reports `(mirror-from-partner unrealizable at
xmin)` on both `gauss_legendre(8)` and `lebedev(17)`. So the identity-vs-mirror
confusion that has bitten twice in this campaign cannot arise for the wrap —
the two maps have disjoint codomains. Record it so a future reviewer does not
spend the budget re-deriving it.

### D.3 Config-blindness ledger (AGENT.md §0.6)

| blindness | applies? | resolution |
|---|---|---|
| flat flux | N/A — no solve here; these are operator-level claims | — |
| 1-group | N/A — no eigenvalue claim anywhere in B3.4c | the trailing probe shape is `(5, 3)` / `(1, 3)`, non-square, so a group↔spatial transpose cannot hide |
| homogeneous | N/A | — |
| **slab is degenerate** | **YES** — a slab has exactly two faces, so "the partner" and "the other face" coincide and a `face_opposite` that returned "any other face" would pass | **C4's 2-D companion is MANDATORY**; `[M]` on `level_symmetric(4)` the wrong-axis partner is shape-legal |
| isotropic-source snapshot | N/A | — |
| **shared draw across faces** | **YES, the campaign's own §16.4 lesson** | every C4/C5 fixture is independently seeded AND asserts the seeds differ |

---

## E. Mutation table

Each row: the production mutation, the test that MUST red, and whether the
mutation is BLIND to some gate that might otherwise be credited. Mutations are
applied by **monkeypatch in-process** — never `git checkout` (the working tree
carries uncommitted state; see `.claude/rules/process-discipline.md`).

| # | Mutation | Must RED | Also reds | BLIND to (named scope) |
|---|---|---|---|---|
| **M1** | `face_name`: swap the suffix table (`{"min": +1, "max": -1}`) | C1 `test_the_bijection_reproduces_the_five_retired_transcriptions`; `test_angular_trace_space.py::test_face_normals_z_derived_from_axis_names` | most of the SN suite (every face's normal flips) | ⚠ **BLIND to C1's round-trip rows** — the round-trip is invariant under relabelling both directions. This is why the `legacy` transcription leg exists. |
| **M2** | `face_normal`: drop the `stem in AXIS_NAMES` check (accept `"minx"`) | C1 negative legs | — | BLIND to every positive row |
| **M3** | `face_opposite`: `return face` (identity) | C2 involution stays GREEN(!); C2 `test_opposite_has_no_fixed_point` reds; C3 `test_the_wrap_consumes_the_OPPOSITE_face` reds; C4 all value rows red | — | ⚠ **BLIND to the involution row** — the identity IS an involution. The no-fixed-point row is what catches it. Record this: involution alone is a weak law. |
| **M4** | `SpatialWrap.domain_face`: `return face` | C3 `test_the_wrap_consumes_the_OPPOSITE_face`; C4 all value rows; B.5 leg 1 | `_assert_wrap_identification` would then compare `Γ₊(face)` vs `Γ₋(face)` → also reds at realization | BLIND to C1/C2 (different object) |
| **M5** | `SpatialWrap.domain_face`: delete the axis-vs-face reconciliation | C3 `test_a_wrap_on_a_foreign_axis_REFUSES` (6 params) | — | BLIND to everything else: `[M]` `PeriodicBoundary(axis="y")` on `xmin` currently produces a *plausible* wrap |
| **M6** | `_assert_wrap_identification`: compare SIZES instead of SETS (`partner_outflow.size == face_inflow.size`) | **NOTHING in the tree today** — closed by the C3 companion below | — | ⚠ **KNOWN-BLIND on every canonical fixture** (`[M]` `\|Γ₊\| == \|Γ₋\|` everywhere ⇒ a size check always passes). **`[M]` VERIFIED CLOSEABLE**, see §E.2. |
| **M7** | `_assert_wrap_identification`: delete it entirely | M6's new negative; C8's missing-partner negative | — | BLIND to every green-path gate — the guard is green by construction on the canonical path (`realizer.py` says so) |
| **M8** | `_reflect_trace`: read `trace.outflow_restriction(face)` instead of the source face (i.e. revert step 4) | C4 `test_the_inflow_is_the_PARTNERS_outflow` (×2), `test_it_is_NOT_the_faces_own_outflow` (×2), C7 `test_B_is_block_structured...`; the 2-D companion | — | ⚠ **BLIND to C5 reciprocity** — `[M]` rel 1.15e-16, exactly reciprocal. Named scope in C5's docstring. |
| **M9** | `_reflect_trace` transpose leg: scatter over `γ₊(face)` instead of `γ₊(partner)` | C5 `test_euclidean_reciprocity_on_a_periodic_slab` (`[M]` rel 4.3e-2 GL8 / 1.5e-2 lebedev17); `test_the_transpose_scatters_over_the_PARTNERS_outflow` | — | BLIND to every forward-only gate |
| **M10** | `_reflect_trace` transpose leg: scatter over `γ₋(partner)` | C5 reciprocity (`[M]` rel 2.3e-1 / 1.0e-1); the support gate | — | — |
| **M11** | Apply the wrap on ONE direction only (xmin reads xmax; xmax echoes itself) | C4's `[face,partner]` parametrisation — the `("xmax","xmin")` param | — | ⚠ **BLIND to a single-face gate.** This is why C4 is parametrised over both directions, not written once for xmin. |
| **M12** | `SpatialWrap.is_adjointable → False` | C6 `test_the_wrap_declares_an_honest_transpose`, `test_every_deck_transformation_is_adjointable` | — | ⚠ **BLIND to `B.is_adjointable` and to every realized-operator predicate** — `[M]` both answered `True` even while the factor said `False`. Named scope in C6. |
| **M13** | `PeriodicWrapOperator.apply → 2.0 * x` | A-3; C4 value rows; B.5 leg 2; A-10b | — | BLIND to C1/C2/C3/C6 |
| **M14** | `PeriodicWrapOperator.apply → x[::-1]` (a permutation) | A-10b (`lebedev(17)`, DISCRIMINATING); C4 on `gauss_legendre` | — | ⚠ **would be BLIND on `level_symmetric` / `product(2,4)`** if those were the fixture — see §D.1. |
| **M15** | `SNMesh.BOUNDARY_OPERATOR_REGISTRY` gains `"periodic"` | A-19 (`test_registry_contains_only_vacuum_and_reflective`, `test_error_lists_supported`) | — | This mutation is the DECISION in 0.1(1). If B3.4c takes it, A-19 migrates and `_law_from_tag` needs a `PeriodicBoundary(axis=AXIS_NAMES[label.axis_index])` arm (explorer §1.3's latent orientation bug) — **and that arm needs its own per-face gate**, or every face gets `axis="x"`. |

### E.1 Mutations that are BLIND everywhere — stated, not hidden

* **`TANGENTIAL_EPS → 0.0`** — blind to the periodic identification on every
  quadrature (§C9), because it moves both half-traces symmetrically. A blind
  mutation NAMED is a scope statement.
* **M3 vs the involution row** — `face_opposite → face` IS an involution, so
  `test_opposite_is_an_involution` stays green. The no-fixed-point row is the
  catcher. Recorded because "involution" reads like a complete characterisation
  and is not.
* **M8 vs C5 reciprocity** — `[M]` rel 1.15e-16. Named in C5's own docstring.
* **M12 vs every realized-operator predicate** — `[M]` both `True` already.

### E.2 `[M]` The M6 closer — VERIFIED to bite, and buildable today

The only input that separates a SET comparison from a SIZE comparison is a
hand-built method space whose `Γ₋` is same-sized and wrong. Measured against
the live guard:

```
F2 — proposed M6 negative: a SAME-SIZED but WRONG inflow set
  honest space: inflow = [4, 5, 6, 7]  outflow = [0, 1, 2, 3]
  honest realize -> OK (guard is green by construction)
  poisoned space: inflow = [0, 1, 2, 3] (same SIZE, wrong SET)
  poisoned realize -> BoundaryError(law='periodic')
      A periodic law on face 'xmin' identifies it with partner 'xmax', but Γ₊(xmax) =
      [4, 5, 6, 7] is not Γ₋(xmin) = [0, 1, 2, 3]. The wrap is realized as the identity ...
F3 — the SIZE-only variant would pass on the poisoned space?
  |inflow| poisoned = 4  |Gamma+(xmax)| = 4  -> a SIZE comparison would PASS: True
```

So the gate is exactly:

```python
def test_the_identification_compares_SETS_not_sizes() -> None:
    r"""The anti-Mode-12 leg of `_assert_wrap_identification`.

    `[M]` `|Γ₊| == |Γ₋|` on every quadrature × face in the tree, so a size
    comparison would be green on every canonical fixture — the error class
    sits inside the shape functional's invariance group. The ONLY input that
    separates the two is a same-sized WRONG set, which no mesh produces and a
    hand-built method space does. `[M]` a size-only guard PASSES on this
    input; the set guard raises.
    """
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    space = face_method_space(quad, face="xmin", faces=("xmin", "xmax"))
    poisoned = replace(space, inflow_indices=np.asarray(space.outflow_indices))
    assert np.asarray(poisoned.inflow_indices).size == \
           np.asarray(space.inflow_indices).size, "the probe must be same-SIZED"
    with pytest.raises(BoundaryError) as exc:
        SNBoundaryRealizer().realize(PeriodicBoundary(axis="x"), poisoned)
    assert exc.value.law == "periodic"


def test_the_honest_space_realizes(self) -> None:
    """The POSITIVE half (vv anti-pattern #11: never test a contract-validator
    only against a broken instance). `[M]` green by construction on the
    canonical path — which is the guard's own stated design."""
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    space = face_method_space(quad, face="xmin", faces=("xmin", "xmax"))
    SNBoundaryRealizer().realize(PeriodicBoundary(axis="x"), space)
```

`[M]` the same guard runs cleanly on the tangential fixture:

```
F4 — product(2,4): does the identification guard run cleanly there?
  xmin: realize OK  |G-|=2 |G+|=2  apply -> (2, 3)
  xmax: realize OK  |G-|=2 |G+|=2  apply -> (2, 3)
```

### E.3 `[M]` The §B.5 re-pose is stateable today, and the endomorphism persists

```
F1 — law.geometry_map.domain_face('xmin') = xmax
     law.geometry_map.domain_face('xmax') = xmin
F5 — realized op type: TensorProductOperator
     realized op.is_adjointable: True
     SpatialWrap('x').is_adjointable: True
     op.apply(ones((|G+|,3))).shape = (4, 3)
     op.apply(ones((N,3))).shape   = (8, 3)   <- still an ENDOMORPHISM (step 4 not done)
```

So `test_b3_domain_narrowing.py[periodic]` (A-11) is still correctly xfailing
against the live tree: leg A passes by accident (`(4,3)` in, `(4,3)` out), leg
B — "does not emit N rows for a full-face input" — still fails at `(8,3)`.
**A-11 flips only when step 4 lands AND the realized wrap validates its own
domain.** If step 4 leaves the leaf accepting a full-face input, A-11 stays red
and B3.4c is not done.

---

## F. `--runxfail` audit

```
$ .venv/bin/python -O -m pytest tests/geometry/test_bc_equivalence_snapshot.py --runxfail -q
...........................................F                             [100%]
=================================== FAILURES ===================================
____ TestPeriodicLebedev17Snapshot.test_delivers_the_partner_faces_outflow _____

    @pytest.mark.xfail(strict=True, reason=_B34C_XFAIL)
    def test_delivers_the_partner_faces_outflow(self) -> None:
        ...
        case = case_by_id(self.case_id)
        snapshot = _load_snapshot(self.case_id)
        space = face_method_space(
            case.build_quadrature(), face=case.face, faces=case.faces,
        )
        actual = _realize(case, space).apply(snapshot["psi_out"])
>       np.testing.assert_array_equal(actual, snapshot["psi_in"])
E       AssertionError:
E       Arrays are not equal
E
E       Mismatched elements: 735 / 735 (100%)
E       First 5 mismatches are at indices:
E        [0, 0, 0]: 1.3969355235725678 (ACTUAL), 0.8927739862960189 (DESIRED)
E        [0, 0, 1]: 1.7590833470217022 (ACTUAL), 0.38813485542612347 (DESIRED)
E        [0, 0, 2]: 0.3349167448547625 (ACTUAL), 0.8763664744929012 (DESIRED)
E        [0, 1, 0]: 1.7428979970872203 (ACTUAL), 1.4820932656512429 (DESIRED)
E        [0, 1, 1]: 1.85048932800639 (ACTUAL), 1.8840142286380854 (DESIRED)
E       Max absolute difference among violations: 1.98492445
E       Max relative difference among violations: 1421.31168824

tests/geometry/test_bc_equivalence_snapshot.py:790: AssertionError
=========================== short test summary info ============================
FAILED tests/geometry/test_bc_equivalence_snapshot.py::TestPeriodicLebedev17Snapshot::test_delivers_the_partner_faces_outflow - AssertionError:
1 failed, 43 passed, 1 warning in 1.97s
```

**CONFIRMED.** The row reds on its OWN `np.testing.assert_array_equal` at
`test_bc_equivalence_snapshot.py:790` — a VALUE failure, 735/735 elements
(100 %) mismatched, max abs diff 1.98. Not a setup error, not a fixture error,
not an incidental exception. The three supporting statements (`case_by_id`,
`_load_snapshot`, `face_method_space`) all completed; the LIVE siblings
(`test_the_two_faces_are_a_periodic_pair`, `test_the_two_probes_discriminate`,
`test_the_identity_body_is_correct_on_the_partner_half_trace`) are inside the
`43 passed`. The Mode-8 class-4 discipline the body claims holds.

The other four strict xfails in the boundary suites, for completeness:

```
$ .venv/bin/python -O -m pytest tests/geometry tests/sn/operators/test_b3_domain_narrowing.py \
    tests/numerics/test_periodic_wrap_operator.py tests/numerics/test_operator_capability_predicates.py -q -rxX
XFAIL tests/geometry/test_bc_equivalence_snapshot.py::TestPeriodicLebedev17Snapshot::test_delivers_the_partner_faces_outflow - B3.4c — periodic is a TWO-FACE coupling ...
XFAIL tests/geometry/test_boundary.py::test_wave0_sum_of_realized_bcs_acts_as_weighted_sum - B3.4 — a NARROWED law cannot be summed with an UN-NARROWED one ...
XFAIL tests/geometry/test_boundary.py::test_operator_sum_of_bcs_acts_as_weighted_sum - B3.4 — a NARROWED law cannot be summed with an UN-NARROWED one ...
XFAIL tests/geometry/test_reemission_closure.py::TestExactlyOneFactorIsNonTrivial::test_exactly_one_of_g_r_is_non_trivial[reflective_a07] - B5 — ReflectiveBoundary(axis, α<1) carries TWO non-trivial factors ...
XFAIL tests/sn/operators/test_b3_domain_narrowing.py::TestEveryLawsDomain::test_realized_law_has_gamma_plus_domain[periodic] - B3.4c — PERIODIC still realizes FULL-FACE ...
633 passed, 4 skipped, 5 xfailed, 5 warnings in 1.67s
```

**B3.4c's todo list is exactly TWO of these five**: the snapshot row and
`test_b3_domain_narrowing.py[periodic]`. The two `test_boundary.py` mixed-law
rows are MISATTRIBUTED (§A-22) and the `reemission_closure` row belongs to B5.

---

## G. Sequencing — the commits, and the gate that must be red before each

Updated after §0.2 — **all seven production steps have now landed on disk; what
remains is entirely the TEST surface**, which is the honest statement of where
B3.4c is.

| step | production | gates it owes | status |
|---|---|---|---|
| **G1** cleanup (1-2) | ✅ on disk | C1, C2, C2b | ⛔ **UNWRITTEN** — `face_name`/`face_opposite`/`FACE_NAMES` have zero coverage; `[M]` all rows green + M1/M3 teeth measured |
| **G2** `domain_face` (3) | ✅ on disk | C3 (incl. the Protocol-completeness leg) | ⛔ **UNWRITTEN**; `[M]` all rows green, M4 reds 25 |
| **G3** realization guard (5) | ✅ on disk | E.2 SET-vs-size pair, C9 `product(2,4)`, C8 re-posed | ⛔ **UNWRITTEN**; `[M]` E.2 verified to bite (probe F) |
| **G4** `is_adjointable` (6) | ✅ on disk | C6 + its scope note | ⛔ **UNWRITTEN**; `[M]` M12 reds exactly 2 rows |
| **G5** ⭐ `_reflect_trace` (4) | ✅ on disk | C4 (1-D **+ the mandatory 2-D**), C5 pair, C7 zero-base | ⛔ **UNWRITTEN**; `[M]` all green, M8 reds 6 |
| **G6** the two live xfail markers | — | §B.5 re-pose of A-10a; `_LAWS["periodic"] → False` in A-11 | ⛔ **OPEN** — `[M]` both still xfail AFTER step 4 |
| **G7** migrations | — | A-1, A-2 (were RED, `[M]` now repaired by the main agent), A-22 (`[M]` repaired) | ✅ partially done |
| **G8** docs (7) + C7 rescoping | partial | falsehoods 1-3 (`numerics/operator.py:2616-2641`, `periodic.py:35-45`), the three block-diagonality claims + their three test docstrings | ⛔ **OPEN** |

**The single most important line in this table: G5's production is IN and G5's
gates do not exist.** Right now the entire cross-face contract rests on two
xfail markers that (§B.2, measured) cannot detect it.

**Retest set** (`.claude/rules/vv-testing.md`, canonical invocation, SERIAL):

```
python -O -m pytest tests/geometry tests/numerics/test_face_layout.py \
  tests/numerics/test_face_name_bijection.py tests/numerics/test_angular_trace_space.py \
  tests/numerics/test_periodic_wrap_operator.py tests/numerics/test_operator_capability_predicates.py \
  tests/sn/operators tests/sn/primitives/test_face_name_crosswalk.py \
  tests/sn/sweep/test_sweep_acyclicity.py tests/transport tests/diffusion -q
```

plus, because step 1 touched `sweep_schedule.py` and `transport/method.py`,
`tests/sn/sweep -q` before the commit. The whole tree at `-m "not slow"` is
≈52 min; the above is the minimum honest set.

---

## G.1 Artefacts

* **This plan** — `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/b34c_verification_plan.md`
* **The dry-run module** (every §C gate, verbatim, runnable) —
  `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/b34c_proposed_gates_dryrun.py`.
  Run from the repo root:
  `.venv/bin/python -O -m pytest scratch/b34c_proposed_gates_dryrun.py -q`.
  It is a SCRATCH artefact: split it into
  `tests/numerics/test_face_name_bijection.py`,
  `tests/geometry/test_boundary_domain_face.py` and a
  `TestPeriodicIsACrossFaceCoupling` class in
  `tests/sn/operators/test_sn_boundary_operator.py`, then delete it.
* **The mutation runner** —
  `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/mutations.py` (+ `_mutplug_*.py`).
  Monkeypatch-only, by design; re-run it after any gate edit.
* **Blast radius** (explorer, complementary) — `scratch/b34c_periodic_blast_radius.md`.

---

## H. Open questions for the plan of record (NOT for me to decide)

1. **Does B3.4c admit `"periodic"` to `SNMesh.BOUNDARY_OPERATOR_REGISTRY`?**
   If NO, step 4 is unreachable-from-a-tag and its only exerciser is the
   injection fixture in §C4 — acceptable, but say so. If YES, `_law_from_tag`
   needs a `PeriodicBoundary` arm or every face gets `axis="x"` (explorer §1.3),
   and A-19 migrates in the same commit.
2. **Does `_reflect_trace(faces=…)` stay exact for a wrap face**, or does it
   raise? The docstring's exactness argument (`boundary.py:388-392`) is a
   splitting claim, and an un-gated one is the ERR-056 shape. §C7 second row.
3. **Should `FaceLabel.face_name` delegate to `face_layout.face_name`?** §C2b.
   Either way the agreement wants the gate.
4. **Where does `test_the_identity_body_is_correct_on_the_partner_half_trace`
   live after §B.5?** Keeping both rows in the same class with near-identical
   bodies invites a future reader to delete "the duplicate" — the one that is
   actually the composition gate. Rename, or move the composition claim out to
   `test_sn_boundary_operator.py` entirely (my recommendation).
5. **Should the REALIZER also refuse a wrap whose partner face is absent**, or
   is the composite's permutation certification the right (and only) home?
   `[M]` realization is silent today. Composite-only is defensible — the law is
   fine, the CONFIGURATION is not — but "realize succeeded" reads as
   "installable" and the sphere case proves it is not. §0.2.1(2).
6. **Does the 2-D cross-axis companion go in this phase?** The defence-in-depth
   here is currently ONE deep: `SpatialWrap.domain_face` refuses a foreign axis,
   so a cross-axis map is unreachable *through the law*. It is NOT caught by
   anything downstream — `[M]` a cross-axis partner is SHAPE-LEGAL
   (`|Γ₋(xmin)| = |Γ₊(ymax)| = 12` on `level_symmetric(4)`) while the index sets
   differ, and `_face_domains`'s permutation certification is on face NAMES, so
   `{xmin→ymax, ymax→xmin, xmax→ymin, ymin→xmax}` passes it. If the axis guard is
   ever relaxed (e.g. for #178's `SymmetryBoundary`, non-opposite gluing), there
   is no second net. **A 2-D `TestPeriodicIsACrossFaceCoupling` companion —
   periodic on the x-pair, reflective on the y-pair, asserting the y-images are
   byte-identical to a control — is the cheap closer and is MANDATORY regardless
   (slab is the degenerate two-face case: `face_opposite` returning "any other
   face" passes every 1-D row).**
