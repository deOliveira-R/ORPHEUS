---
name: issue-257-s10b-chi-mix-production-weighting
description: S10b gate spec — production-weighted convex-average chi_mix over ALL fissile isotopes in compute_macro_xs (replaces first-fissile shortcut). CRITICAL gate-4 finding — ZERO multi-fissile k_eff/flux pins in the suite (no re-baseline needed). Flat-flux honest-scope caveat (Mode-7). U_238.h5 present.
metadata:
  type: project
---

# #257 S10b — production-weighted χ_mix (PRE-IMPL gate spec)

Branch `feature/field-typed-operator-algebra`, HEAD `1d4df67` (S10a
`2b57083` landed + V&V matrix regen). SUT = the χ assignment in
`compute_macro_xs` (`orpheus/data/macro_xs/mixture.py:211-214`). NO
production code written by this agent — spec + skeletons only.

## What S10b changes (the SUT)

TODAY (mixture.py:211-214), a documented "simplification":
```python
chi = np.zeros(NG)
if fissile_indices:
    chi = isotopes[fissile_indices[0]].chi.copy()   # ONLY isotope 0
```
This ignores every fissile isotope past the first. For UO2 (U235+U238,
both flagged fissile, `fissile_indices=[0,1]`) U238's emission spectrum
is dropped entirely.

THE FIX — production-weighted convex average over ALL fissile isotopes:
```
prod_i  = aDen[i] * np.sum(isotopes[i].nubar * sigF[i])   # flat-flux production
w_i     = prod_i / sum(prod_j over fissile_indices)        # Σ w_i = 1
chi_mix = sum(w_i * isotopes[i].chi for i in fissile_indices)
```
All ingredients (`sigF (n_iso,NG)`, `nubar (NG,)`, `aDen`,
`fissile_indices`) are in scope at mixture.py:168-194; `SigP` at :195-197
already sums `nubar*sigF*aDen` per fissile isotope.

## ⭐⭐ Mode-7 HONEST-SCOPE CAVEAT — flat-flux representative weighting (BAKE INTO THE SPEC)

This is NOT the truly-correct material χ. The exact spectral weight is
FLUX-weighted: `w_i ∝ Σ_g νΣf_{i,g} φ_g`. A static data-prep value
CANNOT capture φ. S10b pins the **flat-flux production-weighted convex
average** — i.e. `w_i ∝ aDen_i · Σ_g νΣf_{i,g}` (φ_g ≡ 1). The gate
MUST declare this in every relevant docstring: it verifies the
flat-flux weighting, NOT flux-exactness. (Per the user this is a
documented caveat, NOT a separate issue — do NOT file one, do NOT spec
a flux-exactness leg.) This makes it a flux-shape claim only at the
representative level; eigenvalue-correctness of the resulting χ is
NOT claimed (Mode-7: the simplification is declared, not hidden).

## ⭐⭐ GATE-4 CRITICAL FINDING — ZERO multi-fissile k_eff/flux pins exist

The brief asked me to "enumerate every test with a tight k_eff/flux pin
on a MULTI-fissile mixture … classify which must be re-baselined." The
honest answer after a full audit: **there are NONE.** No re-baseline is
required. Audit trail (HEAD 1d4df67):

- Multi-fissile `compute_macro_xs` callers (the ONLY paths the S10b
  branch is reached on): exactly the three recipes in
  `orpheus/data/macro_xs/recipes.py`:
  - `uo2_fuel` — `fissile_indices=[0,1]` (U235+U238) → MULTI-fissile.
  - `pwr_like_mix` — `fissile_indices=[0,1]` (U235+U238) → MULTI-fissile.
  - `aqueous_uranium` — auto-detect; isotopes `[H_001, O_016, U_235]`
    → only U235 fissile → SINGLE-fissile → χ_mix == U235.chi (gate-3
    byte-identity, NOT a value change).
- `pwr_like_mix` consumers: **zero** (grep `pwr_like_mix` across
  `tests/ orpheus/` returns only its `def`).
- `uo2_fuel` consumers: exactly one — `tests/sn/operators/test_solver_components.py:556`
  inside the `solver_421g` fixture. That fixture feeds ONLY
  `TestPerformanceBaseline::test_profile_421g` (and `test_profile_components`
  uses `solver_2g`, not 421g). `test_profile_421g` is a PROFILING test:
  it prints timings (`print(f"  [421g] ...")`) and asserts NOTHING about
  k_eff or flux. So the one multi-fissile mixture that reaches a solver
  never has its converged value pinned.
- `examples/` demos call recipes but `pyproject.toml` sets
  `testpaths = ["tests"]` → examples are NEVER pytest-collected → no pin.
- `tests/data/test_cross_section_data.py` imports only `_AMU_TO_G` /
  `_number_density` from `recipes` (the number-density helper) — no
  multi-fissile recipe, no χ pin.

⇒ Gate-4 deliverable: **the re-baseline list is EMPTY.** No committed
test reads a converged k_eff/flux off a multi-fissile mixture, so the
χ_mix value change cannot move any green assertion. This is the
inverse of the S10a precursor's L20 scope-blowup: S10a's value-neutral
zeroing touched 377 symbols because non-fissile chi is consumed nowhere
(χ·0); S10b's value-CHANGE touches the converged surface nowhere
because the multi-fissile mixtures are only ever PROFILED, never PINNED.
The verification weight therefore rests entirely on the new
hand-reference gate (gate 2), not on a regression diff.

⚠ Recommend the main agent confirm there are no NON-pytest consumers
(a notebook, a CI smoke script outside `tests/`) that pin a uo2_fuel /
pwr_like k_eff. The grep covered `tests/ orpheus/ demos/ scripts/
examples/`; if a downstream analysis harness lives elsewhere it would
be the sole re-baseline candidate. Treat the EMPTY list as
"none found in the repo tree", confirm before crediting.

## Gates (5)

Test file: `tests/data/test_chi_mix_production_weighting.py`
(`@pytest.mark.foundation` — the convex-average law has NO theory
`:label:`; foundation carries NO `verifies(...)` per [[feedback_vv_tagging]]).
The simplex-preservation interlock (gate 1) REUSES S10a's
`EmissionSpectrum.assert_normalized` verbatim — a correctly-weighted
convex combination of simplices is a simplex, so S10a's guard fires
for free on the S10b output and a weighting bug that breaks the simplex
is caught at construction.

### Gate 1 — convex-simplex intrinsic property (L11, vv #11 interlock)
χ_mix is a probability simplex (Σ=1, χ≥0) for ANY multi-fissile mixture
because a convex combination of simplices is a simplex. Build (HAND,
not via the production expression) several 2+-fissile cases with
DIFFERENT simplex χ_i and varied production weights; assert
`np.isclose(chi_mix.sum(), 1.0, atol=1e-12)` AND `np.all(chi_mix >= 0)`.
Then assert the S10a interlock: `EmissionSpectrum(chi_mix).assert_normalized()`
MUST NOT raise (the same call `Mixture.__post_init__` makes — confirms
the two interlock and the S10b output is a legal `Mixture.chi`). Vary
weights across cases (equal, skewed, near-degenerate) so the simplex
claim is not an accident of a single weight vector.

### Gate 2 — 2-fissile hand-reference (THE TEETH, L11 structurally independent)
A mixture with 2 fissile isotopes whose χ DIFFER and whose production
weights are NON-trivial, so the answer ≠ isotope-0's χ AND ≠ a 50/50
mean. Hand-compute `w_i` and `χ_mix` with an INDEPENDENT formula (NOT
the code's expression) and assert the code matches to FP
(`np.testing.assert_allclose(chi_mix, expected, atol=1e-12, rtol=0)`).
This is the row that distinguishes production-weighting from BOTH the
old first-isotope shortcut (would give χ_0) AND a naive unweighted mean
(would give (χ_0+χ_1)/2). Discriminator design: pick `prod_0 ≠ prod_1`
(e.g. via aDen and/or νΣf so w = [0.25, 0.75] or [0.8, 0.2]) and
`χ_0 ≠ χ_1` so all three candidate answers (χ_0, mean, weighted) are
pairwise distinct — the test mutation-pins the weighting specifically.
- ⭐ STRUCTURAL-INDEPENDENCE of the reference: hand-compute via a
  SEPARATE arithmetic path — e.g. lay the per-isotope production as
  explicit scalars `p0 = aDen0 * (nubar0 * sigF0).sum()`,
  `p1 = aDen1 * (nubar1 * sigF1).sum()`, `chi_mix = (p0*chi0 + p1*chi1)
  / (p0 + p1)` written out term-by-term in the test (NOT a loop that
  re-spells the code's `sum(w_i * chi_i)`). The reference and the SUT
  share only `numpy` arithmetic primitives (below the trusted-library
  line per algebra-of-record); they do NOT share the production
  `compute_macro_xs` expression (above the line).
- ⭐ DRIVING THE SUT: build two minimal SYNTHETIC `Isotope`s and drive
  the REAL `compute_macro_xs`. Single-`sig0` isotopes take the cheap
  paths: `interp_xs_field`/`interp_sig_s` return `field[0].copy()` for
  `n_sig0==1` (interpolation.py:31,64), `_interp_sigT` returns
  `sigT[0,ig]` (sigma_zeros.py:78). So `sigF` passes through VERBATIM
  and `χ_mix` is built from the exact `nubar`/`sigF`/`chi` the test
  controls. CAVEAT: `solve_sigma_zeros` still loops all NG=421 groups
  (cheap for n_sig0=1 but O(NG)); keep the synthetic isotopes at NG=421
  (the `Isotope.NG` constant) with only a handful of nonzero groups so
  `(nubar*sigF).sum()` is hand-traceable. The synthetic fissile isotope
  MUST carry a valid simplex χ (S10a `__post_init__` enforces it on
  construction) — set `chi` to a 421-vector summing to 1 with the mass
  in 1-2 groups. Pass explicit `fissile_indices=[0,1]` to skip
  auto-detect ambiguity.
- ⭐⭐ LOAD-BEARING (PROVEN by live probe at HEAD 1d4df67): every
  synthetic-isotope `compute_macro_xs` drive MUST pass `n_legendre=1`.
  The synthetic isotopes carry one Legendre order (`sigS=[[csr_matrix]]`,
  the S10a `_isotope` shape) and the default `n_legendre=3` indexes
  `iso.sigS[1]` → `IndexError`. The REAL-data smoke leg keeps the
  default 3 (real isotopes carry 3 orders). The skeleton bakes
  `n_legendre=1` into all four synthetic drives.
- ⭐⭐ TEETH-DRY-RUN PROVEN ON REAL HEAD CODE (no stub needed — the
  current first-fissile shortcut IS the wrong formula): with the gate-2
  controlled inputs `aden=[2,1]`, `(nubar,sigf,g)=(2.5,0.10,0)` and
  `(3.0,0.20,5)`, `chi0=[.8@0,.2@1]`, `chi1=[.3@5,.7@6]`, the
  CURRENT `compute_macro_xs` returns `χ_mix == χ_0` (True), `== S10b
  expected` (False), `== unweighted mean` (False). So gate 2 is RED at
  HEAD and goes GREEN exactly when production-weighting lands — and it
  distinguishes the weighted answer from BOTH the shortcut AND the mean.
  `expected = (p0·χ0 + p1·χ1)/(p0+p1)` with `p0=2·2.5·0.10=0.50`,
  `p1=1·3.0·0.20=0.60` → sum 1.0, nonzero groups [0,1,5,6]. The teeth
  are dry-run-verified WITHOUT a stub.
- ⭐⭐ RECOMMENDED ELEGANCE NIT for the implementer (coding-elegance
  Pattern 3 + 5): the χ-mix arithmetic is an unnamed inline sum today.
  Recommend the production change extract a NAMED helper
  `production_weighted_chi(isotopes, sigF, aDen, fissile_indices) -> EmissionSpectrum`
  (the per-isotope `production_rate` is a named reactor-physics
  intermediate `[1/(cm³·s)]` at flat flux; the convex weights `w_i` are
  named). IF the implementer ships that helper, gate 2 can ALSO unit-test
  it directly (faster, no 421-group `compute_macro_xs` drive) AND the
  end-to-end drive becomes a thin integration leg. The helper is NOT
  required for the gates to work (the gates drive `compute_macro_xs`
  end-to-end either way) — flag it as the elegant option, let the
  method-implementer / elegance-enforcer decide.

### Gate 3 — single-fissile byte-identity
A mixture with exactly ONE fissile isotope → `χ_mix == that isotope.chi`
EXACTLY (the weighted average collapses to identity, w=[1]). Bit-identical:
`np.testing.assert_array_equal(mixture.chi, np.asarray(isotope.chi))`.
Proves the change is scoped to multi-fissile and leaves every
single-fissile case untouched. TWO legs:
- (3a) synthetic single-fissile drive of `compute_macro_xs` (mirrors
  gate 2's machinery with one fissile isotope + one non-fissile) →
  `chi == fissile_iso.chi` byte-id.
- (3b) the `fissile_indices=[]` / no-fissile path → `chi == np.zeros(NG)`
  byte-id (the SUT's `chi = np.zeros(NG)` default branch, unchanged by
  S10b). This is the borated_water / zircaloy_clad shape.
Note `np.asarray(...)` on the RHS because `isotope.chi` is now an
`EmissionSpectrum` (S10a); the values must be byte-equal, the subclass
wrapping is irrelevant.

### Gate 4 — k_eff re-baseline enumeration → EMPTY (see GATE-4 FINDING above)
Deliverable is the audit, not a re-baseline. The single-fissile /
synthetic / zero-fissile cases are UNAFFECTED by construction
(byte-identical, gate 3). No multi-fissile converged-value pin exists,
so there is nothing to re-baseline. Document the empty list so the
implementer does NOT go hunting for snapshots to regenerate (a
regenerate-snapshots reflex here would be wrong — there are no
multi-fissile snapshots, and regenerating single-fissile ones would
mask gate-3 if it ever regressed). [[feedback_principled_over_bit_identical]]
still governs IF a downstream multi-fissile pin is later added: it
would be re-baselined as principled-correct (the new χ_mix is more
correct — U238 contributes), NEVER silently.

### Gate 5 — Mode-8 (every assertion -O-safe)
ORPHEUS canonical invocation is `python -O` → bare `assert` STRIPPED.
ALL legs route through `np.testing.assert_allclose` /
`np.testing.assert_array_equal` / `pytest.raises` / a module-level
`_require(cond, msg)` helper (`pytest.fail` — a CALL, fires under -O).
NEVER a bare `assert`. Recommend a TEETH dry-run: stub the production
χ-mix with the WRONG formula (first-isotope shortcut, or unweighted
mean) and confirm gate 2 reddens under `-O` before crediting it; revert
the stub. The single-fissile gate-3 leg also pins the shortcut would-be
behavior is correct THERE (single fissile → both formulas agree), so
gate 2's multi-fissile discriminator is the sole catcher of the shortcut.

## Optional real-data smoke check — U_238.h5 IS PRESENT

Both `orpheus/data/micro_xs/U_235.h5` (49.9 MB) and `U_238.h5` (37.8 MB)
are present. The optional real-data smoke leg is feasible:
- `skipif` guarded on BOTH `.h5` present (mirror S10a's
  `TestRealGendfConstructs` skipif idiom).
- Build a real 2-fissile UO2-like mixture via `uo2_fuel(temp_K=900)` OR
  directly `compute_macro_xs([U235, U238, O16], densities,
  fissile_indices=[0,1])`. Assert: (a) `mixture.chi` is a simplex
  (`assert_normalized` no-raise — the S10a interlock on REAL data); (b)
  `mixture.chi` differs from BOTH `U235.chi` and `U238.chi` (the value
  actually moved — `not np.allclose(mixture.chi, U235.chi)` — proving
  the multi-fissile weighting fired on production data, NOT the
  shortcut). Do NOT pin the exact χ_mix VALUE against a hard-coded
  number on real data (that would be a brittle data-version pin); pin
  the STRUCTURAL facts (simplex + moved-off-both-inputs). This is a
  smoke leg, L11 reference is gate 2's synthetic hand-calc, not this.
- ⚠ `uo2_fuel` prints + runs 421-group sigma-zeros → slow (mark `slow`
  if it exceeds a few seconds; the synthetic gate-2 drive is the fast
  primary). Consider `compute_macro_xs([U235,U238,O16], …)` directly to
  skip the recipe's pyXSteam/density ceremony.

## Regression command (route around baseline reds #250/#232/#212)

Run AFTER production lands, under `-O` (canonical). The S10b touch is
isolated to `compute_macro_xs` χ; the relevant green-anchors are the
data layer + the FissionOperator/keff path (which is BLIND to the χ
VALUE on single-fissile/synthetic fixtures — those use xs_library
`make_mixture`, which builds `Mixture` DIRECTLY and never reaches
`compute_macro_xs`, so they are S10b-inert and MUST stay byte-id):

```
.venv/bin/python -O -m pytest \
  tests/data/test_chi_mix_production_weighting.py \
  tests/data/test_chi_invariant_enforcement.py \
  tests/data/test_emission_spectrum.py \
  tests/data/test_mixture.py \
  tests/data/test_cross_section_data.py \
  tests/sn/operators/test_fission_operator.py \
  tests/sn/operators/test_fission_kernel_crosscheck.py \
  tests/sn/verification/analytical/test_kinf_homogeneous.py \
  tests/sn/regression/test_dd_regression.py \
  -p no:cacheprovider \
  -k "not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical)"
```

- The `-k` excludes the 5 stale SPHERE snapshots (#250, main baseline-red
  w/ #232).
- If the keff suite is widened: `--deselect tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff`
  (#212 `continuous_get` hang).
- These reds are PRE-EXISTING on the branch (NOT S10b) — confirm at
  clean HEAD before crediting any red to S10b.
- ⭐ The DD regression + fission-operator + kinf anchors verify the χ
  VALUE change is INERT on the snapshot/synthetic mixtures (they are all
  xs_library single-fissile region "A" or non-fissile B/C/D, NONE
  multi-fissile, NONE through `compute_macro_xs`). They MUST stay
  byte-identical — if any snapshot moves, the S10b change leaked into a
  path it should not touch (STOP and investigate). The χ_mix VALUE is
  pinned ONLY by gate 2's hand-reference, not by these anchors.

## Anti-recommendations (honored)
- NO flux-exactness leg (flat-flux weighting is the documented
  approximation; pin only the convex-average-with-flat-flux-weights —
  Mode-7 honest scope).
- NO change to single-fissile result (gate 3 byte-identity).
- NO new field type (χ_mix stays a plain (NG,) array that
  `Mixture.__post_init__` coerces to `EmissionSpectrum`).
- NO re-baseline hunt (gate 4 empty — there are no multi-fissile pins).
- NO snapshot regeneration (single-fissile snapshots are S10b-inert;
  regenerating would mask a future gate-3 regression).

## Self-improvement
- NO new vv failure mode. Gates are vv #11 contract-validation
  (gate 1 interlock), L11 hand-reference (gate 2), byte-identity
  (gate 3), Mode-7 honest-scope declaration (the flat-flux caveat —
  already a skill row), Mode-8 (gate 5) — all already in the skill table.
- NO ERR (next free ERR-064 — no real bug; S10b is a documented
  simplification being upgraded, not a caught bug). The OLD first-fissile
  shortcut was a known/documented limitation, not a hidden bug, so it
  does not earn an ERR entry; if the upgrade ever uncovers a latent
  multi-fissile correctness issue in a consumer, THAT earns the ERR.

Extends [[issue-257-s10a-emission-spectrum-verification]] (the S10b
SEAM was pre-declared there — this is its closeout) +
[[feedback_principled_over_bit_identical]] (the value-change-vs-byte-id
governance) + [[feedback_vv_tagging]] (foundation NO verifies) +
[[feedback_regression_tolerance_design]].
