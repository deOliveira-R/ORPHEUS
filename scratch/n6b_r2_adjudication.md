# #340 N6b / R2 — adjudicating the 20 tests that will start warning

`[M]` 2026-08-10, branch `refactor/operator-strategy-layers`, HEAD `0d38d331`.
Population source: `scratch/issue_340_truncated_inner_population.md`.

> Line numbers below are current at `0d38d331` and drift; the durable address
> is the symbol name (`_solve_fwd`, `_TRUNCATED_INNER_CASES`, the frozen
> `_ISO_SPHERE_*` literals). Re-derive if the tree has moved.

**Method.** Every row was MEASURED, not inferred. Each configuration was
reproduced in-process under `$CLAUDE_JOB_DIR/tmp/probe_*.py` and its
`solution.history.record.report()` read — the record already carries the level
tree, the binding criterion, ρ and the projection, so no guard flip and no
monkeypatch of `orpheus/sn/solver.py` was needed. For every row that looked like
it might be pinned by the truncation, a **counterfactual solve at a converging
budget** was run and the test's own assertion re-evaluated against it. No test
file was edited.

---

## 0. Headline

| verdict | rows | which | what it means |
|---|---|---|---|
| **DELIBERATE** | 5 | 14, 17–20 | truncation is load-bearing for the fixture → R2 declaration + narrow `filterwarnings` |
| **INHERITED** | 5 | 7–11 | frozen snapshot captured WITH the truncation; 2 of the 5 **hard-fail** on a budget raise |
| **INCIDENTAL** | 10 | 1–6, 12–13, 15–16 | budget is just short → **no suppression**; list for a follow-up |

⚠ Rows 17, 18 and 20 are counted DELIBERATE but carry an INCIDENTAL half that
must be fixed before the declaration lands (below). Two findings dominate this
document:

- **The four `test_inner_tol_bias_collapses_at_1e_12` rows are only HALF
  deliberate.** Each runs TWO solves. The **loose** leg's truncation is the
  fixture (correct). The **tight** leg — the *reference* the assertion is made
  against — ALSO truncates in three of the four rows, and that is incidental.
  `filterwarnings` is per-TEST, so a bare marker silences both. Fix the tight
  leg FIRST, then declare.
- **`test_cache.py::test_collision_cache_invariance_under_source_iteration` is
  deliberate-by-consequence with a FALSE stated rationale.** Its own
  `len(keff_history) >= 5` floor FAILS at a converging budget (measured: 3
  outers), so the truncation is load-bearing — but the docstring credits those
  ~7 outers to the heterogeneous fixture, and that is measurably wrong.

---

## 1. The adjudication table

ρ and `ran/budget` reproduce the population file. `proj(first)` is
`record.first_failure.projected_iterations()` as the brief asked; `proj(worst)`
is the max over every truncated level in the tree (the first failure is rarely
the worst); `sufficient` is the MEASURED max inner count once nothing truncates.

| # | node-id | verdict | structural reason | budget if INCIDENTAL |
|---|---|---|---|---|
| 1 | `tests/sn/eigenvalue/test_keff_2d.py::TestSIKrylov2DEquivalence::test_eigenvalue_jacobi_gauss_seidel_equivalence` | **INCIDENTAL** | Asserts Jacobi≡G-S to 1e-8; nothing studies truncation. Jacobi truncates **6/6**, G-S 3/6. The docstring's stated mechanism ("the inner SI stops at a slightly different point", i.e. at `inner_tol`) is measurably FALSE — the Jacobi inner never reaches `inner_tol` at all. | proj(first) 1838, proj(worst) 2066 — both LOW; measured sufficient **2601** ⟹ `max_inner=2700` |
| 2 | `tests/sn/eigenvalue/test_keff_curvilinear.py::TestCylinderMultiGroupMultiRegion::test_2g_heterogeneous_fuel_moderator` | **INCIDENTAL** | Asserts only `isfinite`, `k>0`, `0.5<k<3.0`. Budget falls ~60 % short of what `inner_tol=1e-10` costs; 2/3 inners truncate. | proj **788**; sufficient **789** ⟹ `max_inner=800` |
| 3 | `…::TestCylinderMultiGroupMultiRegion::test_2g_heterogeneous_product_different_resolutions` | **INCIDENTAL** | Asserts `\|k(4×8)−k(8×8)\| < 0.05`. BOTH legs truncate 2/3. Raising the budget moves k by 5.5e-11 / 5.8e-11 — nine orders below the bar. | proj **788** (4×8) / **787** (8×8); sufficient 789 / 788 ⟹ `max_inner=800` |
| 4 | `…::TestCylinderMultiGroupMultiRegion::test_multigroup_eigenvector_not_flat` | **INCIDENTAL** | Asserts a fuel/moderator spectrum-ratio difference > 0.01. Same config as row 2. | proj **788**; sufficient **789** ⟹ `max_inner=800` |
| 5 | `…::TestMultiGroupMultiRegionSpherical::test_2g_heterogeneous_converges` | **INCIDENTAL** | Asserts `isfinite`, `0.1<k<3.0`. 3/4 inners truncate. Docstring literally says "must converge to finite keff" — it does not converge. | proj **972**; sufficient **973** ⟹ `max_inner=1000` |
| 6 | `…::TestMultiGroupMultiRegionSpherical::test_multigroup_eigenvector_not_flat` | **INCIDENTAL** | Spectrum-ratio smoke assertion; same config as row 5. | proj **972**; sufficient **973** ⟹ `max_inner=1000` |
| 7 | `tests/sn/regression/test_dd_regression.py::test_dd_regression[2d_2g_LS4_dd_8x4_het_si]` | **INHERITED** | Reproduces the frozen `.npz` **bit-identically today** (`\|Δk\|=0.0`, `max\|Δφ\|=0.0`); the snapshot was generated through the SAME `run_case`/`run_config` path, i.e. with the truncation in place. 7/7 inners truncate. | — |
| 8 | `…::test_dd_regression[cyl_2g_3reg_folded_4x8_dd_n40]` | **INHERITED** | Bit-identical to the snapshot. Raise ⟹ `\|Δφ\|=8.2e-14 ≠ 0` ⟹ `DriftWarning` ⟹ the module's own documented `-W error::DriftWarning` strict PR gate fails. | — |
| 9 | `…::test_dd_regression[slab_2g_3reg_dd_n40]` | **INHERITED** | Bit-identical. ⛔ Raise ⟹ `\|Δk\| = 1.81e-11` vs the correctness gate `SAFETY×keff_tol = 1e-11` ⟹ **hard FAIL**, not merely a re-baseline. | — |
| 10 | `…::test_dd_regression[sphere_2g_3reg_dd_n40]` | **INHERITED** | Bit-identical. Raise ⟹ `\|Δk\|=1.45e-13` passes but `\|Δφ\|=4.1e-13 ≠ 0` ⟹ `DriftWarning`. | — |
| 11 | `…::test_dd_regression[sphere_2g_homogeneous_dd_n20]` | **INHERITED** | Bit-identical. ⛔ Raise ⟹ `\|Δk\| = 2.45e-11` vs the 1e-11 gate ⟹ **hard FAIL**. | — |
| 12 | `tests/sn/solve/test_sn_adjoint_certification.py::TestP13KEquality::test_heterogeneous_slab_k_equality` | **INCIDENTAL** | The warning comes from the shared `_solve_fwd` ONLY (the adjoint leg is fully converged, 0/5). Exactly ONE inner truncates, missing tol by a factor **1.25** (`1.247e-10` vs `1e-10`). | proj **504**; sufficient **505** ⟹ `max_inner=550`. `[M]` k_eff is **bit-identical** at 600 — the fix is numerically free |
| 13 | `…::TestP13Mutations::test_streaming_no_reversal_shifts_k_heterogeneous` | **INCIDENTAL** | Same shared `_solve_fwd`; the assertion is a `>1e-6` mutation shift, so a 1e-10 bias is irrelevant. One edit at the helper fixes rows 12 AND 13. | same as row 12 |
| 14 | `tests/sn/sweep/core/test_cache.py::test_collision_cache_invariance_under_source_iteration` | **DELIBERATE** *(by consequence — see §5: the verdict is a judgment call and the alternative is named)* | 7/7 inners truncate at 50/50. ⛔ **Counterfactual: at `max_inner=425` the power iteration converges in 3 outers and the test's OWN `len(keff_history) >= 5` assertion FAILS.** The truncation is load-bearing — but the docstring's stated reason for the ~7 outers (heterogeneity) is FALSE; they come from the starved inner suppressing the outer increments. | — |
| 15 | `tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py::test_homogeneous_reflective_sphere_iso_unchanged` | **INCIDENTAL** | Looks INHERITED (two frozen in-module literals) and is NOT. `[M]` at `max_inner=350` the solve is fully converged and both literals still pass: `\|Δk\|` 1.87e-11→**5.8e-12** (tol 1e-9, 170× headroom), `max\|Δφ\|` 5.27e-11→**4.55e-11** (tol 1e-9, 22× headroom). The gate's tolerance sits at the FP-tail scale, far above the truncation bias. | proj **313**; sufficient **348** ⟹ `max_inner=350` |
| 16 | `…::test_unclamped_sphere_flux_strictly_positive[2g_R2_S64]` | **INCIDENTAL** | Asserts only `isfinite` and `>0`, but its docstring claims "the **converged** sphere scalar flux" and 4/4 inners truncate. The sibling param `[1g_R2_S16]` does NOT truncate ⟹ any edit must be param-scoped. | proj(first) 2429, proj(**worst**) 2782 ≈ exact; measured sufficient **2790** ⟹ `max_inner=2900`; `[M]` costs only **1.6×** (12.0 s → 19.1 s) |
| 17 | `tests/sn/verification/analytical/test_kinf_homogeneous_tolerance.py::test_inner_tol_bias_collapses_at_1e_12[2eg-cylinder]` | **DELIBERATE (loose) + INCIDENTAL (tight)** | Loose leg (`max_inner=100`, `inner_tol=1e-9`) truncates 4/7 — that IS the studied bias. **The tight leg also truncates 2/5**, and it is the reference the assertion is made against. | tight leg: proj **347**, worst **365** ⟹ `max_inner=420` |
| 18 | `…::test_inner_tol_bias_collapses_at_1e_12[2eg-sphere]` | **DELIBERATE (loose) + INCIDENTAL (tight)** | Loose 4/8 truncated; tight 2/6 truncated. | tight leg: proj **375**, worst **409** ⟹ `max_inner=420` |
| 19 | `…::test_inner_tol_bias_collapses_at_1e_12[4eg-cylinder]` | **DELIBERATE** (clean) | Loose 3/6 truncated; **tight leg is fully converged (0/5)** — the only one of the four that is already clean. This row is purely the intended fixture. | — |
| 20 | `…::test_inner_tol_bias_collapses_at_1e_12[4eg-sphere]` | **DELIBERATE (loose) + INCIDENTAL (tight)** | Loose 4/8 truncated; tight 2/5 truncated. | tight leg: proj **336**, worst **376** ⟹ `max_inner=420` |

**Independent confirmation of the 17–20 split.** The population file records **2**
calls for rows 17, 18, 20 and **1** call for row 19. That is exactly the measured
pattern: two solves warn where the tight leg also truncates, one where it does
not. The call counts corroborate the reading with no extra instrumentation.

---

## 2. DELIBERATE / INHERITED — the prose to paste

### 2.1 Rows 17–20 · `test_inner_tol_bias_collapses_at_1e_12`

⚠ **Sequence matters.** `@pytest.mark.filterwarnings` is per-TEST, not
per-solve, so the marker cannot distinguish the two legs. Land the tight-leg
budget raise FIRST (`max_inner=300 → 420` at
`tests/sn/verification/analytical/test_kinf_homogeneous_tolerance.py:135`),
so the surviving warning is exactly the fixture. `[M]` the raise is safe:
`drift_tight` measures 6.3e-14 … 8.2e-14 against the 1e-11 bar (120–160×
headroom), and `[4eg-cylinder]` already runs its tight leg to convergence at
300 with no ill effect.

**Marker** — function level, on `test_inner_tol_bias_collapses_at_1e_12`
(all four params are deliberate, so the function is the narrowest scope that
covers them without listing params):

```python
@pytest.mark.filterwarnings(
    "ignore::orpheus.numerics.convergence.ConvergenceWarning"
)
@pytest.mark.verifies("matrix-eigenvalue", "fission-matrix", "removal-matrix")
@pytest.mark.parametrize("coord", ["sphere", "cylinder"])
@pytest.mark.parametrize("ng_key", ["2eg", "4eg"])
def test_inner_tol_bias_collapses_at_1e_12(coord, ng_key):
```

⛔ **A second, independent correction is owed on this docstring, found while
adjudicating.** It states *"at default `inner_tol=1e-9` the drift is ~1e-7"*.
`[M]` 2026-08-10 the measured `drift_loose` is **3.9e-11 … 7.5e-11** across the
four rows — **four orders smaller** than the documented figure. The contrast the
test rests on is real but is ~3 orders (7.5e-11 loose vs 8.2e-14 tight), not the
~4+ the prose implies. Fix the number in the same edit; a stale magnitude in a
bias-model docstring is exactly the kind of claim a future session would design a
tolerance against.

**Docstring paragraph to add** (after the existing first paragraph, with the
`~1e-7` in that paragraph corrected to the measured range):

> ⭐ **The LOOSE leg's truncation is the FIXTURE, and it is declared** (#340
> R2). `max_inner=100` at `inner_tol=1e-9` is deliberately starved: the
> under-converged inner IS the bias this test exists to exhibit, so the loose
> solve exits truncated and says so. `[M]` 2026-08-10, 3–4 of its 6–8 inners
> hit the cap at ρ ≈ 0.89–0.93; converging it would need `max_inner ≈ 190–320`
> and would collapse `drift_loose` onto `drift_tight`, deleting the contrast the
> test is built on (`drift_loose` measures 3.9e-11 … 7.5e-11 today against a
> `drift_tight` of 6.3e-14 … 8.2e-14). The suppression above therefore names the
> ONE category this row expects rather than silencing the module.
>
> ⚠ The suppression is per-TEST and this test runs TWO solves, so it covers
> the TIGHT leg as well — which is why the tight leg is budgeted to actually
> converge. `[M]` until 2026-08-10 `max_inner=300` left the tight inner
> truncated in three of the four parametrisations (`2eg-sphere` 2/6,
> `2eg-cylinder` 2/5, `4eg-sphere` 2/5; only `4eg-cylinder` was clean), so the
> *reference* half of a bias measurement was itself biased. `max_inner=420`
> covers the measured worst case (409) with headroom; `[M]` `drift_tight` is
> unmoved at 6.3e-14 … 8.2e-14 against the 1e-11 bar.

### 2.2 Row 14 · `test_collision_cache_invariance_under_source_iteration`

**Marker** — function level (it is a module-level function with a single
`@pytest.mark.l0`):

```python
@pytest.mark.l0
@pytest.mark.filterwarnings(
    "ignore::orpheus.numerics.convergence.ConvergenceWarning"
)
def test_collision_cache_invariance_under_source_iteration() -> None:
```

**Docstring: REPLACE the final paragraph.** The paragraph beginning "A
heterogeneous 2G case is used deliberately…" contains a present-tense-false
claim and must not survive the edit — it credits the ~7 outers to the fixture's
heterogeneity, which is measurably not the cause.

> A heterogeneous 2G case is used deliberately: a homogeneous-reflective slab
> is a degenerate eigenvalue problem (flat flux, `k = νΣ_f/Σ_a` is
> shape-independent), so the fuel|moderator|fuel slab is what gives the flux a
> genuine spatial + spectral shape for the Picard loop to settle.
>
> ⭐ **The truncation is the FIXTURE, and it is declared** (#340 R2).
> `max_inner=50` at `inner_tol=1e-8` is far short of what this slab costs —
> `[M]` 2026-08-10 all 7 inners hit the cap at ρ ≈ 0.958, needing
> `max_inner ≈ 423` to converge — and **the `>= 5` outer floor below depends on
> that starvation.** A starved inner suppresses the outer increments (the #340
> F2 mechanism), which is what stretches this solve to 7 outer steps; `[M]` at
> `max_inner=425` the SAME heterogeneous slab converges in **3** outers and the
> floor assertion FAILS. The invariant under test (`_build_count == 1`) is
> unaffected either way — σ_t is bound once at `SNSolver.__init__` — so the
> truncation costs the gate nothing and buys it the long Picard history it
> asserts against.
>
> ⛔ Until 2026-08-10 this paragraph claimed the ~7 outer steps came from the
> heterogeneity. They do not; they come from the starved inner. The distinction
> matters because it says which knob is load-bearing: raise `max_inner` and this
> test goes red on its own non-degeneracy floor, not on its cardinal assertion.

### 2.3 Rows 7–11 · `test_dd_regression`

**Marker — param level, not function level.** Only 5 of the 13 `CASES`
truncate; a function-level marker would silence the other 8 and any future
case. Wrap the truncating names at the parametrize site so a 6th truncating
case announces itself. `[M]` 2026-08-10 the pattern below was verified end to
end on a throwaway module: the `ids=lambda c: c.name` callable still receives
the underlying `SnapshotCase` (ids render unchanged) and the marker suppresses
the named param ONLY — the unmarked siblings keep warning.

```python
#: The cases whose snapshot was CAPTURED with a truncated inner (#340 R2).
#: Not a suppression of convenience — see the module docstring's
#: "Truncated inners are part of the frozen state" section.
_TRUNCATED_INNER_CASES = frozenset({
    "2d_2g_LS4_dd_8x4_het_si",
    "cyl_2g_3reg_folded_4x8_dd_n40",
    "slab_2g_3reg_dd_n40",
    "sphere_2g_3reg_dd_n40",
    "sphere_2g_homogeneous_dd_n20",
})

@pytest.mark.parametrize(
    "case",
    [
        pytest.param(c, marks=pytest.mark.filterwarnings(
            "ignore::orpheus.numerics.convergence.ConvergenceWarning"))
        if c.name in _TRUNCATED_INNER_CASES else c
        for c in CASES
    ],
    ids=lambda c: c.name,
)
def test_dd_regression(case: SnapshotCase) -> None:
```

**Module-docstring section to add** (after the "Independent corroboration"
paragraph — this is a property of the whole snapshot suite, not of one test
body):

> **Truncated inners are part of the frozen state** (#340 R2). Five cases —
> `2d_2g_LS4_dd_8x4_het_si`, `cyl_2g_3reg_folded_4x8_dd_n40`,
> `slab_2g_3reg_dd_n40`, `sphere_2g_3reg_dd_n40`,
> `sphere_2g_homogeneous_dd_n20` — reach `max_inner` on at least one inner
> before `inner_tol` is met, at ρ ≈ 0.93–0.99. That is not a defect in the
> snapshot: the generator and the test share one `run_case` / `run_config`
> path, so the *truncated* trajectory is deterministic and is exactly what the
> `.npz` froze. `[M]` 2026-08-10 all five reproduce their snapshot
> **bit-identically** (`|Δk| = 0.0`, `max|Δφ| = 0.0`).
>
> ⛔ **Raising `max_inner` re-baselines this suite, and for two rows it is a
> hard failure, not a drift.** `[M]` at `max_inner=2400`:
> `sphere_2g_homogeneous_dd_n20` moves `|Δk| = 2.45e-11` and
> `slab_2g_3reg_dd_n40` moves `|Δk| = 1.81e-11`, both at or past the
> correctness gate `SAFETY × keff_tol = 1e-11`. The other three stay inside the
> correctness gate but every one of the five loses bit-identity, so
> `-W error::DriftWarning` — the strict gate this module documents above —
> would fail on all five. Any budget change here is a snapshot regeneration
> with the independent-corroboration discipline re-run, never a tolerance
> nudge. (`2d_2g_LS4_dd_8x4_het_si` is the stubborn one: ρ ≈ 0.994 leaves it
> still truncated at `max_inner=2400`.)

---

## 3. Do NOT suppress these — the INCIDENTAL follow-up

*(Paste-ready for a GitHub issue body. Ten rows. Every one of these is a true
and useful warning: the solve wanted a converged answer and did not get one.
Silencing any of them is the anti-pattern removed at `3f76d651`.)*

### Group A — one shared helper, numerically free (rows 12, 13)

`tests/sn/solve/test_sn_adjoint_certification.py:92-96`, `_solve_fwd`.
Exactly one inner truncates, and it misses `inner_tol=1e-10` by a factor of
**1.25** (`1.247e-10`). Projected 504; measured sufficient **505**.

- Fix: `max_inner=500 → 550` in `_solve_fwd`, and in `_solve_adj` (line 99-103)
  in lockstep to keep the pair symmetric — the adjoint leg's max used is 469,
  so it is unaffected.
- `[M]` k_eff is **bit-identical** at `max_inner=600` (`|Δk| = 0.0`). No
  downstream assertion moves.
- `test_sphere_k_equality` and `test_infinite_medium_triple_equality` share the
  same helper and do not truncate; the raise is safe for them.

### Group B — curvilinear eigenvalue smoke tests (rows 2–6)

`tests/sn/eigenvalue/test_keff_curvilinear.py`, five call sites all spelling
`max_inner=500, inner_tol=1e-10`.

| config | proj | measured sufficient | suggested |
|---|---|---|---|
| CYL `folded_product(4,8)` (rows 2, 3a, 4) | 788 | 789 | `max_inner=800` |
| CYL `folded_product(8,8)` (row 3b) | 787 | 788 | `max_inner=800` |
| SPH `gauss_legendre(8)` (rows 5, 6) | 972 | 973 | `max_inner=1000` |

⚠ **Consider fixing the tolerance instead of the budget.** These tests leave
`keff_tol`/`flux_tol` at their defaults (1e-7 / 1e-6) while demanding
`inner_tol=1e-10` — the inner is being asked to converge **1000× tighter than
the outer needs**, and that is what costs ~800 iterations. Loosening
`inner_tol` to 1e-8 would converge inside the existing 500 and make the suite
faster; the assertions (`0.5 < k < 3.0`, a spectrum-ratio difference > 0.01, a
0.05 resolution agreement) do not come close to needing 1e-10. `[M]` raising
the budget instead moves k by only 5.5e-11 … 6.7e-11, so either route is safe
for the assertions.

### Group C — W1 iso invariance (row 15)

`tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py:275-279`.
Projected 313; measured sufficient **348** ⟹ `max_inner=300 → 350`.

- `[M]` the two frozen in-module literals SURVIVE the raise with room:
  `|Δk|` 1.87e-11 → 5.80e-12 (tol 1e-9) and `max|Δφ|` 5.27e-11 → 4.55e-11
  (tol 1e-9). This row is NOT inherited — measured, not assumed.
- ⚠ **Cross-check before touching it.** This configuration is the same physical
  problem as the DD case `sphere_2g_homogeneous_dd_n20` (sphere, 20 cells,
  R = 2.0, GL-8, 2G mixture A, reflective, `keff_tol=1e-12`, `flux_tol=1e-10`,
  `max_inner=300`, `inner_tol=1e-10`) and the two have OPPOSITE verdicts — the
  same 2.45e-11 shift is 170× inside the W1 gate's 1e-9 and past the DD gate's
  1e-11. Editing the W1 call site is safe; a "fix the sphere-homogeneous config
  everywhere" sweep is not.

### Group D — the two high-ρ rows (rows 1, 16)

Both were assumed expensive and **measured affordable** — the cost estimate you
would get from the budget ratio is wrong, because converging the inner also
*reduces the outer count* (the #340 F2 mechanism, running the other way).

- **Row 1** — `test_keff_2d.py::…::test_eigenvalue_jacobi_gauss_seidel_equivalence`,
  `max_inner=500, inner_tol=1e-10`. Jacobi truncates 6/6 at ρ = 0.9912
  (projected 1838; `[M]` still 2/6 truncated at 2500). **Measured sufficient
  2601** ⟹ `max_inner=2700`; G-S needs only 1380. `[M]` the Jacobi leg costs
  17.7 s converged vs 20.1 s truncated-at-2500, i.e. total inner work ≈ 2.3× the
  present budget-500 run (~8 s for the pair → ~35 s).
  Independent of the budget decision, **the docstring needs a correction**: it
  states the two schedules agree "to within ~`inner_tol`" because "the inner SI
  stops at a slightly different point", and the Jacobi inner never reaches
  `inner_tol` at all. `[M]` the observed agreement is 6.03e-12 at budget 500 and
  6.54e-13 at 2500 — four orders better than the claimed scale, with the same
  verdict either way.
  ⚠ **Same problem as the INHERITED DD case `2d_2g_LS4_dd_8x4_het_si`** (8×4
  fuel|moderator, LS-4, mixtures A/B), differing only in `max_inner` (500 here,
  300 there) — `[M]` the two converge to the identical limit
  (`1.2119481006243886` vs `1.2119481006243888`). Raising the budget at THIS
  call site is safe and independent; the DD snapshot is a separate frozen
  artefact and must not be dragged along (§2.3).
- **Row 16** — `test_w1_clamp_silent_on_flat.py::test_unclamped_sphere_flux_strictly_positive[2g_R2_S64]`,
  `max_inner=800, inner_tol=1e-11`. 4/4 truncated at ρ = 0.9927; proj(worst)
  2782, **measured sufficient 2790** ⟹ `max_inner=2900`. `[M]` the cost is only
  **1.6×** (12.0 s → 19.1 s), because the outer count drops 4 → 3; `flux.min()`
  is unchanged at 1.326291e-01. This one is cheap enough that raising the budget
  is preferable to loosening `inner_tol` — it also repairs the docstring's
  "the **converged** sphere scalar flux" claim rather than re-wording it. The
  sibling param `[1g_R2_S16]` does NOT truncate, so scope the edit to the S64
  param (or accept the raise on both — the S16 leg's cost is unaffected because
  it already converges inside 800).

---

## 4. Caveat that applies to every projected number

`StoppingCriterion.projected_iterations` fits the observed tail and, as its own
docstring warns, **reads LOW when the budget was cut inside the transient**.
Measured against the budgets that actually sufficed:

| row | ρ | proj(first) | proj(worst) | measured sufficient |
|---|---|---|---|---|
| 12/13 adjoint fwd | 0.9558 | 504 | 504 | **505** |
| 15 W1 iso | 0.9281 | 313 | 347 | **348** |
| 2–4 CYL 4×8 | 0.9729 | 788 | 788 | **789** |
| 5–6 SPH GL8 | 0.9788 | 972 | 972 | **973** |
| 16 S64 | 0.9927 | 2429 | 2782 | **2790** |
| 1 2-D Jacobi | 0.9912 | 1838 | 2066 | **2601** |

So the rule of thumb:

- **proj(worst)** is the number to use, not proj(first) — `first_failure`
  returns the FIRST failing level, which is rarely the worst (on
  `2d_2g_LS4_dd_8x4_het_si` first projects 742 while a later inner projects
  2254). The table in §1 reports both.
- ρ ≲ 0.98 ⟹ proj(worst) is exact to within one iteration; round up.
- ρ ≳ 0.99 ⟹ proj(worst) is a lower bound (row 16: 2782 vs 2790, close; row 1:
  2066 vs 2601, **26 % low**). Add margin or measure.
- ⚠ **Do not cost a budget raise by the budget ratio.** Converging the inner
  removes the F2 increment suppression, so the OUTER count drops — `[M]` row 16
  went 4 outers → 3, and a 3.5× budget raise cost 1.6× wall time.

---

## 5. Ambiguity, stated rather than guessed

**One row is a genuine judgment call: row 14, `test_cache.py`.**

The measurement is unambiguous — a converging budget breaks the test's own
`len(keff_history) >= 5` floor — but which remedy is right is not.

- **Reading A (the verdict above, DELIBERATE):** the truncation is load-bearing
  today, so declare it and move on. Cheapest, and it keeps the row green.
- **Reading B (repair, and arguably more honest):** the floor is measuring the
  wrong quantity. The docstring's actual claim is *"no sweep path is secretly
  rebuilding the cache on every iteration"* — a **sweep**-count claim, not an
  outer-count one. At `max_inner=425` the solve performs ≈ 3 × 425 = **1275**
  sweeps versus 7 × 50 = **350** today, so the cardinal invariant is exercised
  ~3.6× harder at the higher budget; only the outer-count proxy fails. Re-found
  the non-degeneracy floor on sweeps (or on `CollisionCache` call count) and the
  budget can be raised with the gate getting stronger.

The discriminator is a question only the owner can answer: **does the invariant
need many OUTER epochs, or many SWEEPS?** σ_t is bound once at
`SNSolver.__init__` and never rebound inside the solve, which argues the epoch
count is irrelevant and Reading B is correct — but that touches the test's
assertion, not just its budget, so it is out of scope for an R2 declaration
pass. Flagged for the user's call.

**Nothing else is ambiguous.** Every other row was decided on a measurement that
either (a) reproduced the frozen artefact bit-identically, or (b) re-ran the
test's own assertion at a converging budget and reported whether it still holds.

---

## 6. Reproduction

Probes (read-only; no production or test file was modified):

- `$CLAUDE_JOB_DIR/tmp/probe_common.py` — the `record.report()` + `first_failure` dump
- `probe_curvilinear.py` — rows 2–6
- `probe_2d.py`, `probe_2d_hi.py`, `probe_last.py` — row 1
- `probe_dd.py`, `probe_dd_cf.py` — rows 7–11 + the `SAFETY × conv_tol` counterfactual
- `probe_adjoint_cache_w1.py` — rows 12–16
- `probe_counterfactual.py` — the cache `>= 5` floor and the W1 frozen literals
- `probe_sufficient.py` — measured sufficient budgets
- `probe_final.py` — rows 17–20 + the slab DD report
