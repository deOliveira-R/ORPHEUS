# P4.9a — verification plan (PRE-carve)

**Authored** 2026-08-28 by `test-architect`, at `5c4f56d7` (post-P4.3).
**Scope** charter rows 1–3 + the `mm_a_in_coeff` handing. Rows 4–5 (the
`StreamingOperator` 4-arg ctor + `SNMesh` shedding) are **P4.9b, out of scope**.

Every `[M]` below carries the command/probe that produced it. Probes live in
`/private/tmp/claude-501/-Users-rodrigo-git-nuclear-ORPHEUS/292c0918-53a4-478b-b616-9156e051dc68/scratchpad/`
(`probe_deg*.py`, `probe_walk*.py`, `probe_form*.py`, `probe_pin.py`,
`mut.py`, `mut2.py` — all in-process monkeypatch, **no production file was
edited on disk at any point**).

---

## 0. HEADLINE — five findings that change the brief

### F1 ⛔⛔ The brief's named numerical canary is a PROVABLE NON-CATCHER, and so is every frozen snapshot

The done-when says *"the aniso curvilinear canary is bit-identical"*. It will
be — **unconditionally**, because it never executes the changed line.

`[M]` `DiamondDifference.update` and `cell_balance_terms` have **exactly one
production execution route in the whole SN solver**: a CYLINDRICAL mesh whose
`n_phi ≡ 2 (mod 4)`. Full 2-G heterogeneous eigenvalue solves, counting spy on
`DD.update` (`probe_walk2.py`):

| config | `DD.update` calls | `cell_balance_terms` calls |
|---|---:|---:|
| SLAB `gauss_legendre(4)` / `(8)` | **0** | 0 |
| SPHERE `gauss_legendre(4)` / `(8)` / `(9)` | **0** | 0 |
| CYL `folded_product(4, 4)` | **0** | 0 |
| CYL `folded_product(4, 8)` | **0** | 0 |
| CYL `folded_product(4, 6)` | **13 760** (all M-M branch) | 13 760 |
| CYL `folded_product(4, 10)` | **13 536** (all M-M branch) | 13 536 |

The mechanism, `[M]` exhaustively over the even ladder `n_phi = 4…34`
(`probe_deg3.py`): the staggered azimuthal circle places a node at `ω = π/2`
**iff `n_phi ≡ 2 (mod 4)`**, giving exactly one bit-exact `η = 0` ordinate per
μ-level (`deg = n_mu`; `deg = 0` otherwise). `n_phi ≡ 0 (mod 4)` ⟹ `deg = 0`,
16 of 16 values. Sphere is `deg = 0` at every `gauss_legendre(N)` **including
odd N**, where `min |μ| = 0.000e+00` exactly — the μ = 0 ordinate still carries
a downstream face, so it is not degenerate.

⭐ Corroborated independently by the tree's own authored comment at
`tests/sn/sweep/curvilinear/test_apply_matvec_cylinder_invariants.py:70-76`
and by `tests/sn/operators/test_g_adjoint_reciprocity.py:189-210`
(`_make_cyl_product` vs `_make_cyl`). This is not a new discovery — it is a
known fact that the frozen-artifact corpus never adopted.

**Consequently every frozen artifact is blind:**

| frozen artifact | its cylinder quadrature | executes the carve? |
|---|---|---|
| `tests/sn/regression/_generate_snapshots.py:155,157,184,186` | `folded_product(4, 8)` / `(2, 4)` | ⛔ NO |
| `walk_matvec_cyl_2g.npz` via `_make_cyl` (`test_g_adjoint_reciprocity.py:177`) | `folded_product(4, 8)` | ⛔ NO |
| `test_curvilinear_aniso_convergence.py` (the brief's canary) | `n_phi ∈ {8, 16}` | ⛔ NO |
| `TestAffineCarveSweepBaseline[CYL]` (`test_affine_carve_baseline.py:170`) | `n_phi = 2·_N_ORD = 8` | ⛔ NO |

⟹ **The bit-identity set for P4.9a does not exist and must be BUILT.** §3 gives
the cheapest construction (one parametrize row, `[M]` 2.1 ms).

⚠ Filter caveat, stated so the census is re-runnable: my `n_phi\s*=\s*(\d+)`
regex over 874 files reports `n_phi=2 ×8`, but
`test_affine_carve_baseline.py:170` is `n_phi=2 * _N_ORD` — the literal `2` is
a multiplier, not a value. Read every hit's expression before counting it.

### F2 ⛔⛔ There is a THIRD live Morel–Montry spelling, and it is arithmetically DIFFERENT

The brief describes one twin (`diamond.py:229`) of one owner
(`closure.py:1329`). `[M]` `grep -rn "1.0 - tau" orpheus/` returns a third
production site the plan's `transport/`-scoped done-when cannot see:

* **Form A** `(ψ̄ − (1−τ)·ψᵃ)/τ` — `closure.py:1329` (owner) **and**
  `diamond.py:229` (the twin the carve deletes).
* **Form B** `τ⁻¹·ψ̄ − ((1−τ)/τ)·ψᵃ` — `cache.py:377` (`mm_a_in_coeff`) consumed
  at `loss_representation/__init__.py:4348-4352` (the non-degenerate scan fast
  path).

They are algebraically identical and **NOT bitwise equal**. `[M]` on the REAL
τ of `folded_product(4, 6)` (`probe_form.py`, 2400 evaluations): bit-equal
**59.13 %**, max |Δ| `1.776e-15`, max **204 ULP**, median 1 ULP. On uniform
random τ ∈ [¼,¾]: 49.5 % equal, max **37 559 ULP**.

⭐ And the sub-family that hides it: bit-equality is **100 %** exactly when
`τ == 0.5` bitwise (then `1/τ = 2.0` and `(1−τ)/τ = 1.0` are exact), and
54–57 % at `τ = 0.5 ± 1 ULP`. `[M]` on the degenerate ordinates τ is ½ *up to
1 ULP* and exactly ½ only sometimes — `fp(4,6)`: **2 of 4** exact;
`fp(4,10)`, `fp(4,14)`, `fp(6,6)`: **0 of all**. A bit-identity gate validated
on one ordinate of `fp(4,6)` is reading a coin (`vv` #31, #13).

**End-to-end cost of choosing Form B for the degenerate path** `[M]`
(`probe_form3.py`, 2-G heterogeneous cylinder eigenvalue):

| config | `|Δkeff|` | rel | max rel `Δφ` | `keff` bit-equal |
|---|---|---|---|---|
| fp(4,6) nx=8 | 2.776e-17 | 1.86e-16 | 3.78e-16 | **False** |
| fp(4,10) nx=8 | 2.776e-17 | 1.90e-16 | 3.72e-16 | **False** |
| fp(6,6) nx=8 | 0.000e+00 | 0 | 3.75e-16 | True |
| fp(4,6) nx=16 | 2.776e-17 | 1.87e-16 | 4.98e-16 | **False** |

⟹ **RULING the plan should adopt: the carve routes the degenerate branch
through FORM A.** It is the owner's spelling, it costs nothing, and it is the
only choice under which "bit-identical" is a claim rather than a hope. Form B
is numerically harmless (1–2 ULP, inside `iteration_count × ULP`) but breaks
`array_equal` on 3 of 4 configs — for zero architectural gain.

⚠ And the honest scope note the charter needs: **P4.9a does not achieve "the
M-M relation has exactly ONE spelling in the tree".** It achieves
`grep -c "1.0 - tau" orpheus/transport/ == 0`. Form B survives at
`cache.py:377`. The two forms partition the ordinate set (A on degenerate, B
on non-degenerate) so they never disagree on one input — a Pattern-2 twin, not
a numerical inconsistency. Say so in the phase's close-out, or the next session
reads the done-when as stronger than it is.

### F3 ⛔ Row 3 DESTROYS `LinearDiscontinuous`'s curvilinear refusal guard, and its only witness

`[M]` `linear_discontinuous.py:138-153` — `_require_slab` keys the #158
curvilinear refusal on **`upstream_state.angular_upstream is not None`**, the
exact field row 3 removes. Its sole witness is
`tests/transport/spatial/test_linear_discontinuous.py:236
::TestLDGuards::test_curvilinear_visit_raises`, which **constructs
`UpstreamState(angular_upstream=...)` directly** — after row 3 that test is
unwritable.

This is `vv` #17's displaced-guard clause and #28's temporal-twin clause at
once. **Row 3 must re-key the guard and re-write its witness in the SAME
step** (`plan-authoring` §6c).

`[M]` the available replacement signal is a VALUE signal already on the visit
(`probe` in §5 of this file): `visit.streaming_terms.delta_A_over_w` is
**exactly `0.0`** on CARTESIAN (4/4 cells) and non-zero on SPHERICAL
(`62.9 … 9.3`) and CYLINDRICAL (`4.29`). Preferred per `vv` #17: a guard keyed
on a value signal is reachable by calling the scheme directly, so its witness
needs no mesh.

### F4 ⭐ The twin IS gated (13 rows) — but by end-to-end gates, never by an equivalence gate

The memo says *"nothing gates them against each other"*. True and narrower than
it reads: nothing compares the two SPELLINGS, but corrupting the twin does red
things. `[M]` in-process mutation `arm=twin`
(`CellResult.outgoing_angular_state *= 1.05` inside `DD.update`):

**3 unit rows** (`tests/sn/sweep/core/test_diamond.py`, 541-test slice, 74 s):
`TestBitIdenticalCurvilinear::test_spherical_outward_bit_identical`,
`::test_spherical_inward_bit_identical`,
`TestCylindricalDegenerate::test_degenerate_cell_synthetic`.

**10 end-to-end rows** (108-test `n_phi=6` slice, 19 s):
`test_apply_matvec_cylinder_invariants.py::test_cylinder_three_way_standoff[6-{2,4}-{10,20,40}]`
(6), `test_psi_half_positivity.py::test_the_production_value_path_keeps_psi_half_POSITIVE[6]`,
`::test_a_zero_seed_drives_psi_half_NEGATIVE_by_the_amplification[6]`,
`test_sweep_inverse_identity.py::TestSweepInverseIdentity::test_forward_of_inverse_is_identity_on_a_random_composite[cyl_folded]`,
`test_loss_transpose_solve.py::test_g3_full_field_solve_reciprocity[cyl_degenerate]`.

So the carve has real exposure and a real (if narrow) net. **All 13 must stay
green through the carve** — they are the acceptance set, not the canary.

### F5 ⛔ Retiring `cell_balance_terms` breaks 4 DECLARED `implements` edges and orphans an equation's only catchers

`[M]` `mcp__nexus__provenance_chain`:

* `cell_balance_terms` declares `implements` on **`dd-curvilinear-scalar`,
  `dd-cylindrical-degenerate`, `dd-slab-scalar`,
  `streaming-action-cell-balance`** (the rest are `inferred` name matches).
* `cell_balance_for_streaming` already declares 3 of those 4 — **`dd-slab-scalar`
  is NOT among them.** Retiring the helper orphans that edge unless migrated.
* `[M]` `math:equation:dd-mm-closure-constants` has **zero declared**
  implementers (all `inferred via: closure`) and exactly **three** claiming
  tests — `test_diamond.py::TestBitIdenticalCurvilinear::{outward,inward}_bit_identical`
  and `TestCylindricalDegenerate::test_degenerate_cell_synthetic`. **All three
  are the twin's catchers.** If the carve retires or re-scopes them without
  replacement, the equation drops from `verified` to `implemented`.
* `[M]` `DiamondDifference.update` itself declares `implements
  dd-cylindrical-degenerate`. After the carve `update` closes only the SPATIAL
  axis, so that claim becomes **false** and must move to the closure.

⚠ No `@pytest.mark.catches("ERR-NNN")` sits on any affected file
(`[M]` grep over the 4 helper files + `tests/transport/spatial/`), so marker
migration is confined to `verifies`.

---

## 1. Deliverable 1 — the twin-equivalence rewires

### 1.0 The memo's list of five is INCOMPLETE — the real anchor is `TestResidual`

`[M]` in-process mutation `arm=cbt` (`cell_balance_terms.denom *= 1.05`), over
the four helper-naming files (140 tests, 0.79 s), **unclipped** counts:

| file | reds |
|---|---:|
| `tests/sn/sweep/core/test_diamond.py` | **33** (25 `TestResidual` + 5 `…_bit_identical` + `test_degenerate_cell_synthetic` + the 2 `…_matches_cell_balance`) |
| `tests/sn/sweep/core/test_cache.py` | **19** |
| `tests/sn/sweep/core/test_cell_balance_for_streaming.py` | **3** |
| `tests/sn/sweep/core/test_ordinate_scan.py` | **0** |
| **total** | **55** |

Two corrections to the brief fall straight out:

* **`test_ordinate_scan.py` reds ZERO.** Its `_visits_to_ab` (`:64`) uses the
  helper to build the `(a, b)` it then FEEDS to `ordinate_scan`. Both sides
  move together — `vv` #22's shared-input blindness. It is a **fixture
  generator, not a reference**, so the rewire is mechanical and costs no claim.
* **`TestResidual` is the load-bearing twin-equivalence family and is not on
  the memo's list.** `test_residual_zero_multi_group_heterogeneous[*]` and
  `test_residual_bit_identity_at_zero_source[*]` red because `update` (via
  `terms`) and `residual` (via `for_streaming`) disagree — i.e. they ARE the
  cross-helper consistency gate, more so than the two named `..._matches_cell_balance`
  rows.

### 1.1 Per-gate prescription

Legend — **PROMOTE**: claim class strengthens (cross-impl → literal anchor).
**PRESERVE**: claim class unchanged. **DEMOTE**: weakens; docstring must say so.

| # | gate | what it still tests after the retirement | prescription | class |
|---|---|---|---|---|
| R1 | `test_cell_balance_for_streaming.py::test_n_mask_1_matches_scalar_curvilinear` (`:160`) | **Nothing** — its entire subject is the two helpers | Rewrite as a **hand-written literal pin** of `cell_balance_for_streaming` on the module's existing synthetic packet (τ=0.75, α=(0.1,0.15), `A_down=1.5`, `Σ_t=[1.2,0.8]`, `ψˢ=[0.5,0.7]`, `ψᵃ=[0.4,0.6]`). Compute the four expected floats OFFLINE by hand/SymPy and paste them; assert `array_equal`. ⛔ Do NOT compute the reference by re-calling the helper. | PROMOTE |
| R2 | `…::test_n_mask_1_matches_scalar_slab` (`:210`) | Nothing (same) | Same treatment on the slab packet. Keep the two `assert_array_equal(angular_*, zeros)` legs — they are independent structural pins that survive. | PROMOTE |
| R3 | `test_cache.py::test_cache_populator_matches_cell_balance_terms` (`:474`) | **Everything it tested.** `[M]` `affine_scan_coefficients` (`diamond.py:640-693`) computes `denom` from its own expression and does **not** call either helper — so cache-populator vs balance-helper remains two independent implementations | Re-point the reference to `cell_balance_for_streaming` (build `angular_denom_term = dA_w·c_out`, `angular_numer_upstream = dA_w·c_in·ψᵃ` from `visit`/closure). Keep `rtol=1e-14`. **Rename** `…_matches_cell_balance_for_streaming`. | PRESERVE |
| R4 | `test_diamond.py::TestResidual::test_curvilinear_residual_matches_cell_balance` (`:1396`) | ⛔ Becomes `f(x) == f(x)`: `residual` ALREADY routes through `cell_balance_for_streaming` (`diamond.py:304`), so re-pointing the reference makes both sides one call | **Hand-arithmetic reference in the test body**: `denom = 2·|μ|·A_down + dA_w·c_out + Σ_t·V`; `numer = |μ|·A_total·ψˢ + dA_w·c_in·ψᵃ`; `ref = denom·cell_avg − (source + numer)`. Written in the SAME association order (bit-identity survives — `[M]` both helpers already associate `(streaming + angular) + collision`). Docstring: *"the reference is now hand-written arithmetic, NOT a second production helper; this gate pins `residual`'s rearrangement against the balance equation as written on the theory page, and can no longer detect a divergence between two helpers because there is only one."* | DEMOTE (declared) |
| R5 | `…::test_cylindrical_degenerate_residual_matches_cell_balance` (`:1440`) | Same as R4 | Same treatment on `_cylinder_degenerate_visit_inputs()`. This row is the **degenerate** arm — keep it, it is one of the few unit-tier degenerate gates. | DEMOTE (declared) |
| R6 | `test_ordinate_scan.py::_visits_to_ab` (`:64`) + its 2 consumers | Its subject is `ordinate_scan`, untouched | Mechanical swap to `cell_balance_for_streaming` at `n_mask=1`. Add a docstring line: *"the helper here GENERATES the `(a,b)` the SUT consumes — it is a fixture, not a reference; `[M]` a mutation of the helper reds 0 of these rows."* | PRESERVE (already blind — now stated) |
| R7 | `TestResidual` — `[M]` **25 rows** red under M3: `test_residual_zero_multi_group_heterogeneous` ×15, `test_residual_zero_at_solved_cell_avg` ×5, `test_residual_bit_identity_at_zero_source` ×5 | Round-trip `residual(update(...)) ≈ 0` — still a real claim about the **solve/apply pair**, now through ONE helper | **Leave the bodies untouched**; add to the class docstring: *"⚠ before P4.9a these rows also cross-checked `cell_balance_terms` against `cell_balance_for_streaming` (`[M]` 55 reds under a `terms.denom` mutation). After the retirement both directions read one helper, so what survives is the solve↔apply rearrangement, not helper equivalence."* | DEMOTE (declared) |
| R8 | `test_diamond.py::TestBitIdenticalCurvilinear::{outward,inward}_bit_identical`, `TestCylindricalDegenerate::test_degenerate_cell_synthetic` | These consume `DD.update`'s **angular** output, which row 1 deletes | See §1.2 — marker migration, not a rewire. |

### 1.2 The external pins that already anchor these values — hunted, and one MOVES

Per `coding-standards` ("hunt for an EXTERNAL hand-written pin before
concluding you lost coverage — and check it MOVES under the old value"):

* ✅ **`test_cell_visit_c_stamp.py`** (3 foundation tests) pins the
  `c_in`/`c_out` **stamp** against an inline surrogate independent of both
  helpers. It survives row 3 only if `CellVisit.c_in/c_out` survive — see
  Open ruling Q2.
* ✅ **`tests/sn/sweep/core/_c_surrogate.py::c_from_constants`** is the shared
  hand-transcribed `(c_in, c_out)` source. Independent of the helpers; keep.
* ⛔ **No hand-written pin exists for `denom`/`numer_upstream` themselves.**
  R1/R2's literal pins are therefore a genuine new anchor, not a formality —
  this is the coverage the retirement would otherwise destroy.
* ✅ **`TestAffineCarveSweepBaseline[SPH|CYL]`** is a real external bit-identity
  pin and **it moves** — `[M]` §5's `mm_ulp` arm reds exactly those 2 rows.

### 1.3 Marker migration

| marker | on | after |
|---|---|---|
| `verifies("dd-curvilinear-scalar")` | R4, R5 | keep (still the balance equation) |
| `verifies("dd-mm-closure-constants")` | the 3 rows of R8 | ⛔ **must move to a closure-side gate.** These are `dd-mm-closure-constants`' ONLY claiming tests (`[M]` nexus). Write the replacement in `tests/sn/sweep/curvilinear/` against `closure.tau/c_in/c_out_per_ordinate` **in the same commit**. |
| `verifies("dd-cylindrical-degenerate")` | `test_degenerate_cell_synthetic` | keep — but the test's subject moves from `DD.update`'s angular output to the walk's. Re-point the assertion, not the marker. |
| `verifies("blelloch-1990-eq-1-5")` (module `pytestmark`) | `test_ordinate_scan.py` | unchanged |
| `sentinel` | `test_spherical_outward_bit_identical` | keep; re-verify it still reddens after the rewire (`sn_sentinel_harness` discipline) |
| `implements` (doc/registry) | `cell_balance_terms` → `dd-slab-scalar` | ⛔ **migrate to `cell_balance_for_streaming`** — it is the one declared edge the survivor does not already carry |
| `implements` | `DiamondDifference.update` → `dd-cylindrical-degenerate`, and its docstring's M-M claims | re-scope to SPATIAL-only; the angular claim moves to `MorelMontryAngularSweep` |

---

## 2. Deliverable 2 — the `is`-identity gate (the §6c witness)

### 2.1 What it must assert, and why the obvious spelling is vacuous

The charter asks for *"the degenerate-ordinate walk branch and the vectorized
branch reach the SAME closure object"*. If the implementation simply reads
`self.mesh.pole_angular_closure` in both branches, an `id()` comparison is
tautological. The gate earns its name by asserting the thing that is FALSE
today:

> **the degenerate branch reaches a `PoleAngularClosureBase` instance at all,
> and it is the mesh's own.**

`[M]` today it reaches `DiamondDifference` — a `DiscretizationScheme` — and
never a closure. `[M]` baseline for the §6c witness (re-verified at
`5c4f56d7`): test files naming `precompute_psi_state` = **8**, naming
`outgoing_angular_state` = **3**, **intersection EMPTY** (both counts are
positive controls for the negative).

⛔ **The gate's buildability CONSTRAINS the design.** It is only writable if the
degenerate branch calls a **closure METHOD**. If the implementation instead
inlines `geom.tau_inv` / `geom.mm_a_in_coeff` (Form B), no closure object is
reached, the gate cannot be written, and F2's bit-identity is lost as well.
Two independent reasons for the same ruling.

### 2.2 Spec

**File** `tests/sn/sweep/curvilinear/test_angular_closure_is_single_object.py`
(new). Sits beside the other closure gates; `tests/sn/sweep/core/` is the
*spatial* scheme's home and this gate is about the angular one.

**Markers** `@pytest.mark.foundation` (a software-invariant gate: object
identity, no theory `:label:`), module-level `pytestmark`. Per
`feedback_vv_tagging`, foundation carries NO `verifies()`.

**Fixture** — the ONLY family that reaches the branch:

```python
_QUAD = Quadrature.folded_product(n_mu=4, n_phi=6)   # n_phi ≡ 2 (mod 4)
# cylinder, nx=8, ng=2, edges 0.01..2.0, bc_left=reflective, bc_right=vacuum,
# heterogeneous sigma_t + fixed-seed random per-ordinate source
```

`[M]` this config has `deg = 4` degenerate ordinates of `N = 12`, and one
`sweep_once` costs **2.1 ms** and enters `DD.update` **32 times**, every one
through the M-M branch (`probe_pin.py`).

**Observation surface** — least-invasive, and it observes the ROUTE rather than
the value (`vv` #26): monkeypatch the *bound class's* angular-advance method
with a recorder that appends `self` and delegates. Patch on
`type(sn.pole_angular_closure)`, not on the instance, so a branch that
resolves the closure by any route is still caught.

**Legs** (all four; the gate is not done with one):

1. **REACHED** — after one `sweep_once`, `len(seen) > 0`.
   ⛔ *This is the leg that is RED today* — record that in the docstring.
2. **IDENTITY** — `all(obj is sn.pole_angular_closure for obj in seen)`.
3. **BOTH BRANCHES** — run a second `sweep_once` on `folded_product(4, 8)`
   (`deg = 0`, `[M]` §F1) and assert the recorder also fires there, from the
   non-degenerate path, on the same instance. Without this leg the gate proves
   only that *one* branch reached a closure.
4. **NON-VACUITY / negative** — build a SECOND `SNMesh` over the same
   `(quad, coord)` and assert **none** of its closure calls appear in the first
   mesh's recorder (i.e. the recorder discriminates instances). Guards against
   an `is`-assertion that would pass against any object.

**Docstring text to commit verbatim (the §6c red-today record):**

```
⛔ RED BEFORE P4.9a, by construction — and that is this gate's whole point.
Before the un-weld the degenerate cylindrical-axis branch of
``_OneDimScanWalk._run`` closed the angular axis by calling
``DiamondDifference.update``, whose Morel-Montry block was an inline
Pattern-2 twin of ``sn/angular/closure.py`` — not a closure OBJECT at all,
so leg 1 (``len(seen) > 0``) could not pass however the gate were written.
[M] 2026-08-28, at 5c4f56d7: 8 test files named ``precompute_psi_state``,
3 named ``outgoing_angular_state``, and the intersection was EMPTY — nothing
in the tree compared the two spellings.  [M] the branch is reached ONLY on a
cylindrical mesh with ``n_phi ≡ 2 (mod 4)`` (one bit-exact ``η = 0`` ordinate
per μ-level); at ``n_phi ≡ 0 (mod 4)`` ``DiamondDifference.update`` is called
ZERO times in a full eigenvalue solve, which is why every frozen snapshot in
the tree was blind to the twin.
```

**Companion structural leg** (row 3's half, same file):
`CellResult` / `UpstreamState` / `CellVisit` no longer carry the angular
members — assert by `dataclasses.fields`, not `hasattr` on an instance
(a defaulted field would still answer `hasattr` after removal via `getattr`
fallbacks; see `coding-standards`' string-form clause).

---

## 3. Deliverable 3 — the bit-identity set

### 3.1 What must be pinned

The per-cell path runs on cylindrical degenerate ordinates. Pre- and post-carve
arithmetic is the same **only under the F2 Form-A ruling**. Under that ruling
the expectation is **bit-identity**, asserted with `np.array_equal` (`vv`
bit-identity criterion — an implementation move, not a math move).

### 3.2 ⭐ The cheapest construction: one parametrize row on the existing harness

`tests/sn/sweep/core/test_affine_carve_baseline.py` is already the right
instrument — a `--capture-baseline` snapshot harness on a raw `sweep_once`
with **heterogeneous σ_t + fixed-seed random per-ordinate source**, i.e. the
regime that ACTIVATES the curvilinear redistribution term (`vv` §H2), which is
exactly what a flat/vacuum gate cannot see. `[M]` its `_GEOMS_1D = ("SLB",
"SPH", "CYL")` with `_N_ORD = 4` ⟹ `n_phi = 8` ⟹ blind.

**Prescription: add a fourth member `"CYL_DEG"` using
`Quadrature.folded_product(n_mu=4, n_phi=6)`** (do NOT change `_N_ORD`, which
would silently re-baseline the existing CYL row).

`[M]` measured cost and activation (`probe_pin.py`):

| row | N | degenerate ords | `sweep_once` | `DD.update` calls |
|---|---:|---:|---:|---:|
| `CYL` (fp(4,8), existing) | 16 | 0 | 10.1 ms | **0** |
| `CYL_DEG` (fp(4,6), NEW) | 12 | **4** | **2.1 ms** | **32**, all M-M |

Two snapshots (`sweep_angular_CYL_DEG`, `sweep_scalar_CYL_DEG`), ~2 ms, and the
carve acquires a real numerical anchor. This is the single highest-value item
in the plan.

### 3.3 The rest of the set

| # | config | what it pins | cost `[M]` | status today |
|---|---|---|---|---|
| B1 | `TestAffineCarveSweepBaseline[CYL_DEG]` — fp(4,6), nx=8, ng=2, het σ_t, random Q | raw `sweep_once` angular+scalar, bit-identical | ~2 ms | **NEW (§3.2)** |
| B2 | The 13 twin-catchers of §F4 | end-to-end behaviour through the carve | 19 s (108-test slice) + 74 s (541-test slice) | exists; must stay green |
| B3 | `test_g_adjoint_reciprocity.py` / `test_loss_transpose_solve.py` `[cyl_degenerate]` rows (`_make_cyl_product`, fp(2,6)) | the transpose/adjoint arm on degenerate ordinates | in B2's 19 s | exists |
| B4 | A 2-G heterogeneous cylindrical **eigenvalue** solve at fp(4,6), nx=8 | `keff` to `array_equal` | `[M]` **0.62 s** | ⛔ **ABSENT** — recommend one row in `tests/sn/solve/` |
| B5 | `walk_matvec_cyl_2g.npz` | the matvec/adjoint frozen baseline | existing | ⛔ **BLIND** (fp(4,8)) — recommend a sibling `walk_matvec_cyl_deg_2g` built from `_make_cyl_product`, which already exists and is already used by two suites |
| C1 | SLAB `gauss_legendre(8)` + SPHERE `gauss_legendre(8)` | **control**: must be bit-identical *and* must NOT execute `DD.update` | in B2 | exists |

⚠ C1 is not decoration. `[M]` `DD.update` is called 0 times on slab and sphere,
so their bit-identity carries **zero** information about the carve — the
control exists to keep a future reader from citing them as evidence.

---

## 4. Deliverable 4 — the mutation plan

**Harness discipline (crash-safety).** Every arm below is an **in-process
monkeypatch installed by a pytest plugin** (`-p mut`), so **no production file
is ever modified on disk** and a `SIGTERM` at the harness timeout cannot leave
the tree mutated. This is strictly stronger than copy-aside + `diff -q`
(`process-discipline`'s clause), because there is nothing to restore. The
plugin used for the pre-measurements is in the scratchpad and can be committed
under `tests/_mutation/` if the team wants it re-runnable.

⚠ Run every arm with `--continue-on-collection-errors` and count `^ERROR`
separately from `^FAILED` — a mutation that makes production raise kills
COLLECTION and pytest then reports `FAILED = 0` (`vv` #17, and my own L47:
6 of 13 arms once read as clean zeros this way).

⚠ Scope each arm to the files that can redden. `[M]` `tests/sn/sweep/core` +
`tests/transport/spatial` = 541 tests / **74 s**; the `n_phi=6` slice = 108
tests / **19 s**; the four helper files alone = 140 tests / **0.8 s**.

### 4.1 The battery

| arm | mutation | scope | expected | `[M]` pre-carve |
|---|---|---|---|---|
| **M0 control** | `cell_balance_for_streaming` → `denom *= 1.05` | core+transport | many reds | **37 reds / 43 s** ✅ harness alive |
| **M1** | owner `closure.py:1329` M-M → `*1.05` | core + `n_phi=6` slice | ⛔ **must red the degenerate rows AFTER the carve; `[M]` today it does NOT reach them** | run post-carve; this is the phase's proof |
| **M1b** | owner M-M, per-arm: perturb ONLY the `(1−τ)` factor | same | reds | post-carve |
| **M1c** | owner M-M, orientation flip `τ → 1−τ` | same | reds | post-carve. ⚠ my L47: the three committed τ properties (`[0,1]`, `[¼,¾]`, `τ_m+τ_{M−1−m}=1`) are ALL invariant under this flip — do not expect the τ gates to catch it; the catcher is a signed law |
| **M2 (pre-carve baseline)** | `DD.update` → `outgoing_angular_state *= 1.05` | core+transport, `n_phi=6` slice | reds | **13 reds** (3 + 10) — §F4 |
| **M3 (pre-carve baseline)** | `cell_balance_terms` → `denom *= 1.05` | 4 helper files | reds | **55 reds** (33/19/3/0) — §1.0 |
| **M4** | delete the degenerate branch's angular write (`psi_angle[:, i] = …`) | `n_phi=6` slice | reds | post-carve |
| **M5** | route the degenerate branch through **Form B** | B1 + B4 | ⛔ must red `array_equal` and NOT red `allclose` | `[M]` `Δkeff = 2.78e-17`, `keff` bit-equal **False** on 3 of 4 configs (§F2) |
| **M6** | LD guard: make `_require_slab` a no-op | `tests/transport/spatial/` | must red the re-keyed witness | §F3 — proves the re-key kept teeth |
| **M7** | `mm_a_in_coeff := tau_inv - 1.0` | `tests/sn/sweep/core` | reds | **2 reds / 56 s** — §5 |
| **M8 anti-twin** | re-introduce the M-M relation inside `diamond.py` | n/a | ⛔ **must not be spellable**: assert by AST/grep that `transport/` contains no `tau` symbol | the done-when's own grep, as a gate |

### 4.2 The proof obligation the phase owes

**M1 is the phase's whole point.** Post-carve, mutating `closure.py:1329` must
red at least the 3 degenerate unit rows AND the `[6-*]` end-to-end rows —
demonstrating that the per-cell path now pins the single spelling. If M1 reds
only the non-degenerate rows, the degenerate branch did not actually route
through the owner and the carve is cosmetic.

⚠ **Read the red SET, not the count** (`vv` #17's identity clause): M1's red
set post-carve must be a **superset** of M2's pre-carve red set. If it is
merely the same size but disjoint, something else changed.

---

## 5. Deliverable 5 — the `mm_a_in_coeff` handing gate

### 5.1 What is being moved and why the naive gate is tautological

`cache.py:373-377` reads `tau = closure.tau_per_ordinate` and derives
`tau_inv = 1.0/tau`, `mm_a_in_coeff = (1.0 - tau)/tau`. Moving those two lines
INTO the closure, over the same `_tau_per_ordinate_cache` array, is
**bit-identical by construction** — same input, same ops, same order. A
before/after comparison of the moved code against itself proves nothing.

`[M]` the current derivation is exactly the stated expression, 4/4 configs
(`probe_coeff.py`): `tau_inv == 1.0/tau` and `mm_a_in_coeff == (1.0-tau)/tau`
both `array_equal` True on fp(4,6), fp(4,8), fp(2,6), fp(6,8).

### 5.2 ⭐ The realistic mutation, and it is a 1–2 ULP one

The thing that will actually go wrong is an implementer writing the *cleaner*
algebraic form inside the closure:

```python
mm_a_in_coeff = tau_inv - 1.0        # algebraically identical, NOT bit-identical
```

`[M]` `mm_a_in_coeff != tau_inv - 1.0` at **1–2 ULP** on all four configs.

### 5.3 The pin — it EXISTS, it is thin, and it MOVES

`[M]` mutation arm `mm_ulp` (`mm_a_in_coeff := tau_inv - 1.0`) over
`tests/sn/sweep/core` (463 tests, 56 s) reds **exactly 2**:

```
tests/sn/sweep/core/test_affine_carve_baseline.py::TestAffineCarveSweepBaseline::test_sweep_angular_and_scalar_unmoved[SPH]
tests/sn/sweep/core/test_affine_carve_baseline.py::TestAffineCarveSweepBaseline::test_sweep_angular_and_scalar_unmoved[CYL]
```

Baseline for the same slice: 463 passed. So the pin is real (it moves under the
old value — `coding-standards`' licence requirement) and it is the **only** one:
two rows for a `(N,)` field consumed on every curvilinear scan.

### 5.4 Prescription

* **Cheapest pin: extend the existing harness, do not write a new one.** §3.2's
  `CYL_DEG` row adds a third and fourth row to the same gate at ~2 ms.
* **Add one direct field gate** in `tests/sn/sweep/core/test_cache.py`
  (`@pytest.mark.foundation`), asserting `array_equal` for the whole `(N,)`
  closure-algebra block against the closure's own τ:

  ```
  assert np.array_equal(geom.tau_inv,       1.0 / closure.tau_per_ordinate)
  assert np.array_equal(geom.mm_a_in_coeff, (1.0 - closure.tau_per_ordinate)
                                            / closure.tau_per_ordinate)
  assert np.array_equal(geom.c_in,  closure.c_in_per_ordinate)
  assert np.array_equal(geom.c_out, closure.c_out_per_ordinate)
  ```

  ⚠ Docstring must state what this gate **cannot** see: once the closure owns
  the derivation, the right-hand sides are the closure's own expression, so
  this becomes a *spelling* pin (it catches `tau_inv - 1.0`, `vv`-measured at
  1–2 ULP) and **not** an independent-value pin. Name M7 as its proof.
* **Negative leg**: assert the four fields on a SLAB mesh are the neutral
  element (`tau_inv == 1.0`, `mm_a_in_coeff == 0.0`, `c_in == c_out == 0.0`),
  `array_equal`. Without it the gate has no input that could fail structurally.
* ⛔ **Do NOT use `allclose`.** `[M]` the realistic defect is 1–2 ULP; any
  tolerance ≥ `1e-15` makes the gate a non-catcher.

---

## 6. Deliverable 6 — gate order

**Step 0 — baselines, before any edit.**
`[M]` already taken at `5c4f56d7`: `tests/sn/sweep/core` + `tests/transport/spatial`
= **541 passed, 1 skipped, 4 xfailed / 75.5 s**; the `n_phi=6` slice =
**108 passed / 25.8 s**; `tests/sn/sweep/core` alone = **463 passed / 56.7 s**.
Re-take on the branch tip before trusting them.

**Step 1 — LAND THE ANCHOR FIRST (§3.2), still on the pre-carve tree.**
Add `CYL_DEG` to `_GEOMS_1D` and capture its two snapshots. This is a
`--capture-baseline` run against **unmodified production**, so the snapshot is
inherited from verified code (`snapshot_migration_when_production_goes_bare`
rule 4). ⛔ If this lands *after* the carve, the snapshot inherits the NEW code
and pins nothing.
*Gate:* the new row passes; the existing 3 rows are byte-unchanged.

**Step 2 — commit the two pre-carve mutation baselines (M2, M3) into the plan
record.** Numbers, not prose: 13 and 55, with their red sets. They are the
denominators every later claim is read against.

**Step 3 — land the `is`-identity gate RED** (§2), `xfail(strict=True,
reason=…)` with the §2.2 docstring. Per `feedback_vv_tagging` and `vv` Mode 8
class 4: strict, so it must flip; and structure the body so **exactly one
statement can fail and it is leg 1**.

**Step 4 — row 1 + row 2 (the carve proper): DD sheds the M-M block;
`cell_balance_terms` retires.**
Land with R1–R7 (§1.1) in the SAME commit — §6b: `cell_balance_terms`'s test
call sites and its production call site are one call-site set, and R1/R2 have
no subject the moment the helper dies.
*Gates:* B1 bit-identical; the 13 rows of §F4 green; 541-slice green;
`grep -c "1.0 - tau" orpheus/transport/` == 0.

**Step 5 — row 3 (protocol sheds the angular members) + the LD re-key (§F3).**
Same commit as the re-keyed `test_curvilinear_visit_raises` — §6c. The
`is`-identity gate's `xfail` is removed here (or at step 4, whichever step
makes the walk call the closure).
*Gates:* M6 reds; `test_discretization_scheme_protocol.py:179` migrated;
`tests/transport/spatial` green.

**Step 6 — the `mm_a_in_coeff` handing (§5), with its gate.**
*Gates:* M7 reds; the 2 existing `TestAffineCarveSweepBaseline` rows + the new
`CYL_DEG` rows all `array_equal`.

**Step 7 — the post-carve battery.** M1, M1b, M1c, M4, M5, M8 + the M0 control.
Record the table, not a boolean (`vv` #17).

**Step 8 — provenance + docs (§F5).** Migrate the `dd-slab-scalar` declared
edge; re-scope `DiamondDifference.update`'s implements claim; write the
`dd-mm-closure-constants` replacement gate. Then `mcp__nexus__dead_references`
(the only instrument that reads the docstring surface), `sphinx -W`, `pyright`.

---

## 7. Open rulings — BLOCKING

**Q1 ⛔ (blocking, §F2). Which arithmetic form does the degenerate branch use
after the carve?** Recommendation: **Form A** (the owner's
`(ψ̄ − (1−τ)ψᵃ)/τ`), because it is the only choice under which the phase's
"bit-identical" done-when is a claim rather than a coincidence, and because the
`is`-identity gate is unbuildable without a closure METHOD call. Cost of Form
B, if chosen anyway: `Δkeff ≈ 2.8e-17` (1–2 ULP) and `array_equal` fails on
3 of 4 configs, so B1/B4/B5 must then be `assert_array_almost_equal_nulp`
with a stated nulp and a re-baseline justification.

**Q2 (blocking for row 3 scope). Do `CellVisit.c_in` / `c_out` survive?**
The charter sheds `CellVisit.{tau, c_in, c_out}`. `tau` must go — it is the
angular closure's own weight and the twin's only input. But `c_in`/`c_out`
enter the **denominator and numerator**, which stay in the spatial balance;
`residual` already consumes them via `angular_denom_term`/`angular_numer_upstream`.
Two options: (a) shed all three and have the walk pass the two assembled
angular contributions to `update` (mirrors `residual`, mirrors
`cell_contribution`'s return type — recommended); (b) shed `tau` only. Option
(a) changes `update`'s signature ⟹ its own §6b call-site set
(`[M]` `.update(` appears in 14 test files; ~5 are DD-relevant). Option (b)
leaves `c_in`/`c_out` on an L2 protocol under angular names, which is the
smell the phase exists to remove.

**Q3. What replaces `dd-mm-closure-constants`' three claiming tests?**
(§F5.) They are its only ones. A closure-side gate on
`tau/c_in/c_out_per_ordinate` against `_c_surrogate.c_from_constants` is the
natural home and needs no mesh.

**Q4. Should `walk_matvec_cyl_deg_2g` be added (B5)?** `_make_cyl_product`
(fp(2,6)) already exists and is already consumed by two suites, so the
generator change is one `CASES` row. It closes the last blind frozen artifact.

**Q5. Where does the anti-twin gate (M8) live?** Recommendation: extend the
existing layer/AST gate family rather than a new file — a grep-based
`@pytest.mark.foundation` row asserting `orpheus/transport/` contains no `tau`
identifier and no `1.0 - tau` spelling. It is the done-when, made permanent.

---

## 8. What this plan does NOT cover

* Rows 4–5 (`StreamingOperator` ctor, `SNMesh` shedding) — P4.9b.
  `[M]` the memo's sizing (~150 ctor sites, ~60 mesh reads) is unverified here.
* `angular_adjoint` — untouched by P4.9a (§5b/O-3's step). The two spy/count
  gates (`test_psi_half_coupling.py:2633`,
  `test_streaming_cell_transpose_relocation.py:161`) survive.
* The §6 correctness item (curvilinear LD moment metric) — separate, #158.
* Whether `D @ T_ang @ T_spatial` is the right factorization — that is ⛔ P0
  and this phase only makes it measurable.
