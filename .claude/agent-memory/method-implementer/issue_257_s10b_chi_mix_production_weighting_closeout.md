---
name: issue-257-s10b-chi-mix-production-weighting-closeout
description: S10b closeout — χ_mix in compute_macro_xs becomes a production-weighted convex average over ALL fissile isotopes (replaces first-fissile shortcut). NEW module-level helper production_weighted_chi; convex avg expressed as a weights·spectra contraction (weights @ fissile_spectra) NOT sum(generator) — the latter cost +1 net-new pyright. 8/8 gates green incl. real-UO2 smoke. Net-new pyright 0, 0 type:ignore. NOT committed.
metadata:
  type: project
---

# #257 S10b — production-weighted χ_mix (closeout)

Branch `feature/field-typed-operator-algebra`, on HEAD `1d4df67` (S10a
`2b57083` + peierls #264 fix `0aabde4` already committed). Built on the
test-architect gate spec
`issue_257_s10b_chi_mix_production_weighting.md`. **NOT committed**
(per brief). Working tree: `orpheus/data/macro_xs/mixture.py` modified,
`tests/data/test_chi_mix_production_weighting.py` un-skipped (was the TA
skeleton, untracked).

## Behavioral change

`compute_macro_xs` χ assignment (was mixture.py:211-214, the documented
first-fissile shortcut `chi = isotopes[fissile_indices[0]].chi.copy()`)
→ production-weighted convex average over ALL fissile isotopes:
`χ_mix = Σ_i w_i χ_i`, `w_i = a_i Σ_g ν̄_{i,g}σ_{f,i,g} / Σ_j (…)`
(flat-flux representative weighting, φ_g≡1). For UO2 (U235+U238 both
fissile) U238's emission spectrum now contributes instead of being
dropped.

## Deliverables

- **Production**: NEW module-level helper `production_weighted_chi(
  isotopes, sigF, aDen, fissile_indices) -> np.ndarray` (coding-elegance
  Pattern 3 named reactor-physics intermediate `production` `[1/(cm³·s)]`
  at flat flux + named convex weights `w_i`; Pattern 5 reusable
  primitive). Returns a PLAIN ndarray — `Mixture.__post_init__` coerces
  to validated `EmissionSpectrum` (coercion belongs to the container,
  not the data-prep helper, per brief). The χ assignment now calls the
  helper. Guard on `total>0` (degenerate ν̄≡0 → uniform weights, never
  divide-by-zero).
- **Gate skeleton un-skipped**: all 8 tests in
  `tests/data/test_chi_mix_production_weighting.py` made green (was
  8 skipped). Gates: (1) convex-simplex intrinsic + S10a interlock,
  (2) 2-fissile hand-reference TEETH, (3a/3b) single-fissile + no-fissile
  byte-identity, real-UO2 smoke.

## ⭐ THE ONE DEVIATION (load-bearing) — convex avg as a contraction, NOT sum(generator)

The brief's literal helper body ended
`return sum(w * isotopes[i].chi for w, i in zip(weights, fissile_indices))`.
That `sum(generator)` idiom (the SAME one the surrounding production
`SigP = sum(...)` / `Sig2 = sum(...)` already uses) cost **+1 net-new
pyright** (`reportReturnType`: builtin `sum` is typed returning `int`,
so pyright infers the helper returns `int` not `ndarray`). Confirmed by
isolation: mixture.py at HEAD = 3 pyright errors; with the
`sum(generator)` helper = 4; brief requires net-new 0.

**Fix (more elegant, not a workaround):** express the convex average as
the weights·spectra contraction it mathematically IS —
```python
fissile_spectra = np.array([np.asarray(isotopes[i].chi) for i in fissile_indices])
return weights @ fissile_spectra   # χ_mix = Σ_i w_i χ_i, types as ndarray
```
The stacked per-isotope simplex matrix `fissile_spectra` (n_fissile, NG)
is a named intermediate, and `weights @ spectra` reads like `Σ_i w_i χ_i`
(coding-elegance master standard: code reads like the math). This is
strictly better than the brief's `sum(generator)`: typed-clean AND more
domain-legible. mixture.py back to 3 errors (== HEAD baseline). NO
`# type: ignore` added. The pre-existing `SigP/Sig2 = sum(...)` lines
were NOT touched (out of S10b scope — they keep their pre-existing
int|ndarray pyright noise; #226 territory).

## Verification (all pasted in the report)

1. Skeleton 8/8 green under `-O` (incl. real-UO2 smoke — U_235.h5 +
   U_238.h5 present; multi-fissile weighting fired on real data,
   χ_mix moved off both U235.chi and U238.chi).
2. **Teeth proven** by mutation: reverting the production line to the
   legacy first-fissile shortcut → gate 2
   (`test_two_fissile_matches_hand_weighted_average`) reddens AND the
   real-UO2 smoke leg reddens (χ_mix == χ_U235), while gate 1 (simplex
   interlock) and gate 3 (single-fissile byte-id) STAY green — exactly
   the TA-predicted discriminator behavior (gate 2 is the SOLE shortcut
   catcher; single-fissile both formulas agree; χ_0 is itself a simplex
   so the interlock can't see the shortcut). Re-confirmed with the final
   contraction form. Restored after.
3. **TA full-regression command**: 115 passed, 2 xfailed, 0 failed
   (route around #250 SPHERE via `-k`). DD/fission/kinf anchors
   byte-identical. PROVEN scoped two ways: (a) grep — NO regression
   anchor file references `compute_macro_xs` (they route through
   `xs_library.make_mixture` / build `Mixture` directly → S10b-inert);
   (b) the identical DD DriftWarnings (314242 ULP scalar_flux, 132 ULP
   k_eff) appear byte-for-byte at HEAD WITHOUT the S10b change (stash
   test) — pre-existing branch FP-noise, all within tol, NOT S10b.
4. **pyright net-new = 0**: full project 2353 errors, 19 warnings
   (== baseline 2353). mixture.py 3 errors (== HEAD baseline). 0 new
   `# type: ignore`.

## Gate-4 re-baseline list — confirmed EMPTY

No multi-fissile k_eff/flux pin exists in the suite (TA audit). The χ
VALUE change is pinned ONLY by gate 2's synthetic hand-reference + the
real-UO2 smoke; no snapshot regenerated, no value assertion moved. (Per
brief: if a green value assertion HAD moved, STOP — none did.)

## Posture (what was deliberately NOT done)

- NO Sphinx `:label:` / theory stub (the helper docstring carries the
  math + the Mode-7 honest-scope caveat in-module; the convex-average
  law is a software invariant pinned by `@pytest.mark.foundation`, no
  theory-page claim — same posture as S5/S6 type-surface carves). NO
  algebra-of-record Branch-1 SymPy manifest (this is a data-prep
  weighting formula + numpy-primitive delegation; verification = the
  L11 hand-reference + the simplex interlock, not a Branch-1 derivation).
- NO new ERR (the first-fissile shortcut was a DOCUMENTED limitation
  being upgraded, not a hidden caught bug — TA self-improvement note;
  ERR-064 stays reserved).
- NO archivist DISPATCH this stage (no rich-narrative claim minted; the
  honest-scope caveat lives in the helper docstring). If a future stage
  mints a theory label for the mixture-spectrum approximation, THAT
  earns the archivist pass.
- NO flux-exactness leg (flat-flux weighting is the declared Mode-7
  approximation, NOT a separate issue per the user/TA).

## LESSON (skill-relevant)

When a brief hands a literal `return sum(generator-of-arrays)` body, the
builtin `sum` types as `int` under pyright and costs +1 net-new
`reportReturnType` — the surrounding production may already carry that
noise (`SigP/Sig2 = sum(...)`) but a NEW symbol must not add to it under
a pyright-net-new-0 gate. The principled fix is NOT `# type: ignore`
nor `np.sum` (wrong axis semantics) — express the reduction as the
linear-algebra contraction it IS (`weights @ stacked_spectra` for a
convex average), which both types cleanly AND reads like the math
(`Σ_i w_i χ_i`). A convex average is a matrix-vector product; spell it
that way. (coding-elegance master standard + Pattern 3: the stacked
spectra matrix is the named intermediate.)

Extends the S10a→S10b seam (S10a `2b57083` minted `EmissionSpectrum` +
source-chi enforcement; S10b is the consumer that produces a
multi-fissile χ_mix the S10a guard validates for free — the gate-1
interlock).
