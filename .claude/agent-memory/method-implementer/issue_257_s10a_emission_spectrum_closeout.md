---
name: issue-257-s10a-emission-spectrum-closeout
description: S10a — EmissionSpectrum(ndarray) self-validating value-object + Mixture/Isotope __post_init__ simplex/null χ guard. ADDITIVE + behavior-NEUTRAL (no χ value change). Precursor blast radius was 3× the brief's scope (9 extra non-fissile direct-construction sites + 1 production-bearing-SigP=0-SigF subtlety). 0 net-new pyright, byte-id DD regression GREEN. REFINED 2026-06-21 (review outcome): guard re-keyed SigF→PRODUCTION (νΣf>0), single-sourced via enforce_emission_spectrum helper; SigF/SigP seam RETIRED (billiard hack reverted); promoted ERR-063 peierls source-index test.
metadata:
  type: project
---

# #257 S10a — EmissionSpectrum newtype + container χ-invariant guard (CLOSEOUT)

Branch `feature/field-typed-operator-algebra`, on HEAD `c6e21c0` (matches
the TA gate-spec probed HEAD). NOT committed. ADDITIVE + behavior-NEUTRAL
precursor: NO χ numerical value changes (that is the deferred S10b
multi-fissile production-weighting). Closes the type+guard half;
`compute_macro_xs` (mixture.py:191-194) still uses the first-fissile-isotope
χ verbatim — the documented S10b seam.

## Deliverables (manifest)

1. **NEW leaf** `orpheus/data/emission_spectrum.py` — `EmissionSpectrum(np.ndarray)`
   subclass. `__new__` = `np.asarray(values).view(cls)`; `__array_finalize__`
   no-op (no carried state). `assert_normalized()` (Σχ≈1 with
   `np.isclose(atol=1e-12, rtol=0)` AND `np.all(self>=0)`, the two clauses
   INDEPENDENT — `[1.2,-0.2]` reds on the negativity clause despite Σ=1);
   `is_emitting` property `bool(np.any(self>0))` (Python bool, NOT np.bool_);
   `assert_null()` (`np.all(self==0)`, EXACT zero — constructed zeros).
   Imports ONLY numpy (cycle-safe — both micro_xs.isotope and macro_xs.mixture
   import THIS). Module + class docstring states: SOURCE-LEVEL validated
   value-object carrying the probability-simplex law; explicitly NOT the
   #263-deferred per-cell `SpectrumField` (no flux-algebra membership); tol
   rationale (real GENDF residual ≈2.2e-16, atol 1e-12 = ~4 orders headroom).
2. **`Isotope.__post_init__`** (`micro_xs/isotope.py`) — coerce
   `self.chi=EmissionSpectrum(self.chi)`, then `if is_fissile:
   assert_normalized() else: assert_null()`. (`is_fissile` already existed,
   keyed on `np.any(self.sigF>0)`.)
3. **`Mixture.is_fissile` property + `__post_init__`** (`macro_xs/mixture.py`)
   — NEW `is_fissile = bool(np.any(self.SigF>0))` (NOTE: I wrapped in `bool()`,
   so it's a Python bool and pyright-clean — Isotope's twin LACKS the wrap and
   carries a PRE-EXISTING `reportReturnType` red at isotope.py:71, untouched).
   Same `__post_init__` pattern as Isotope.

4. **Foundation gates (TA-authored skeletons, made green, ZERO edits to them):**
   `tests/data/test_emission_spectrum.py` (28 legs incl. real-GENDF) +
   `tests/data/test_chi_invariant_enforcement.py`. Both `@pytest.mark.foundation`,
   NO `verifies()` (simplex law has no theory `:label:` equation — it's a
   software invariant). Run `-O`: **28 passed** (the 1 PytestConfigWarning is
   the expected `-O` notice; ALL substantive legs route through
   `_require`/`np.testing`/`pytest.raises`, Mode-8-safe).
5. **Sphinx stub** `docs/theory/cross_section_data.rst` — NEW section
   "The Fission Emission Spectrum as a Validated Value-Object" with
   `:label: emission-spectrum-simplex` (the simplex/null law equation),
   `:mod:`/`:class:`/`:meth:` cross-refs, and a `.. todo:: Archivist
   expansion needed (S10a)` carrying the rich-narrative brief (5 points incl.
   the S10b seam + the SigF-vs-SigP subtlety). Build clean (see below).

## Precursor edits (behavior-NEUTRAL — landed in this same change)

The TA gate-spec's L20 precursor finding was CORRECT but INCOMPLETE. It
named (a) the 9 non-fissile xs_library region dicts B/C/D and (b) 2
`test_sweep_cache.py` fixtures. Those + the strict guard surfaced a WIDER
blast radius: **9 ADDITIONAL non-fissile DIRECT-construction sites** the
factory-path audit could not see (they call `make_mixture(...)`/`Mixture(...)`
DIRECTLY with inline `chi=[1.0,...]` placeholders, NOT via
`get_mixture("B"/"C"/"D")` — so the library zeroing does NOT fix them
transitively).

- `xs_library.py` (per brief): `_B/C/D_1G` `chi=np.array([1.0])`→`np.zeros(1)`;
  `_B/C/D_2G` `[1,0]`→`np.zeros(2)`; `_B/C/D_4G` `[.6,.35,.05,0]`→`np.zeros(4)`.
  Region A (fissile) UNTOUCHED. Behavior-neutral: B/C/D have `sig_f≡0` ⇒
  `νΣf≡0` ⇒ `χ·(νΣf·φ).sum≡0`.
- `test_sweep_cache.py:287/:598` (per brief) → `np.zeros(1)`.
- **EXTRA (found via a deterministic AST scanner, NOT in the brief):**
  - `tests/moc/test_verification.py` `_make_pure_absorber_1g`/`_make_pure_scatterer_1g`
    (non-fissile, `chi=[1.0]`→`np.zeros(1)`). The two MoC reds that first
    surfaced the wider radius (`test_pure_scatterer_equilibrium_single_sweep`,
    `test_scatter_only_source`). [`_make_fission_only_*` + `_make_n2n_*` are
    FISSILE — chi=[1,0]/[1.0] is a valid simplex, LEFT untouched.]
  - `tests/cross_method/cases.py` `_make_unit_sigma_t_one_group_mixture` —
    the `c<1.0` non-fissile branch shared a `chi=np.array([1.0])` with the
    fissile `c>=1.0` branch; split so the non-fissile branch builds
    `np.zeros(1)`, the fissile branch keeps `[1.0]`. (Was a COLLECTION error
    via module-level fixtures in `test_eigenvalue.py`/`test_polymorphism.py`.)
  - `tests/sn/operators/test_scattering_kernel_crosscheck.py` `_mix`,
    `test_scattering_operator.py` `test_n2n_doubling_factor` +
    `test_residual_zero_when_p0_diagonal_only_no_n2n` (non-fissile, `chi=[1,0]`→
    `np.zeros(2)`).
  - `tests/sn/primitives/test_axis_native_construction.py` (COLLECTION error,
    module-level `_MATERIALS`) + `test_axis_primitive.py` `_one_group_mixture`
    (`chi=[1.0]`→`np.zeros(1)`).
  - `tests/sn/solve/test_d3_admission.py` pure-absorber 2g (`chi=[1,0]`→`np.zeros(2)`).
  - `tests/sn/sweep/core/test_sweep_schedule_nd.py` + `test_unified_sweep_dispatch.py`
    (2 sites) pure-scatterer (`chi=[1.0]`→`np.zeros(1)`).

## ⭐ THE ONE TRUE BLAST-RADIUS FINDING (SigF vs SigP semantic gap)

`tests/derivations/test_trajectory_resolvent_billiard.py::_mixture_from_xs`
builds a **MULTIPLYING medium**: `SigP = nu_sigma_f > 0` (e.g. 0.1) but
`SigF = 0` (the helper takes `nu_sigma_f` directly, never separating Σf and ν).
The billiard reference solver (`billiard.py:1031-1032`) reads `mix.SigP` and
`mix.chi` to build the fission source — it NEVER reads `mix.SigF` (nor
`absorption_xs`, grep-confirmed across the whole `trajectory_resolvent`
module). So:
- `is_fissile` (keyed on `SigF>0`, per the brief's pinned definition)
  returns **False** → the guard demands `assert_null()` → `chi=0`.
- But `chi` IS consumed (× `SigP>0`) → zeroing it would CHANGE the multi-group
  billiard result. NOT behavior-neutral. The pure-precursor fix (zero chi) is
  WRONG here.

**Resolution (behavior-preserving):** the material genuinely emits, so its
TRUE classification is fissile. Set `SigF=nu_sf_arr.copy()` in the helper —
this flips `is_fissile`→True so the valid-simplex default `chi=[1,0,...]`
passes `assert_normalized`, and it is INERT for the test because the billiard
path reads only `SigT`/`SigS`/`SigP`/`chi` (never `SigF`, never
`absorption_xs`; `SigT` is passed directly, not recomputed). Documented
in-code as a stand-in (true `Σf=SigP/ν`, ν untracked). The 12 billiard tests
pass.

This is a GENUINE design seam: **`is_fissile = np.any(SigF>0)` mis-classifies
a production-bearing material whose `SigF` field happens to be zero.** For
REAL data (`compute_macro_xs`) it never bites — `SigP=nubar*sigF` so
`SigP>0 ⟺ SigF>0`. It only surfaces in synthetic test helpers that set `SigP`
directly. NOT changing the pinned `is_fissile` definition (out of S10a scope);
flagged for S10b/future consideration (could key `is_fissile` on `SigP` OR
require helpers to set `SigF` consistently). Logged for the user.

## Verification (all four, verbatim)

1. **TA skeletons (`-O`):** `28 passed, 1 warning in 0.41s` (the warning is
   the `-O` PytestConfigWarning — informational).
2. **Byte-identity DD regression (`-O`):** `13 passed, 13 warnings in 45.25s`
   — the 13 DriftWarnings are WITHIN-tol + PRE-EXISTING (TA-confirmed at HEAD,
   NOT an S10a effect). NO snapshot moved ⇒ the precursor zeroing of inert
   non-fissile χ is PROVEN behavior-neutral (a moved snapshot would mean χ was
   NOT inert; none moved).
3. **TA full regression command (`-O`, route around #250/#232):**
   `107 passed, 2 xfailed, 13 warnings in 60.05s` (2 xfailed = pre-existing
   operator-subset, TA-confirmed). PLUS a codebase-wide sweep (all edited-site
   suites): `321 passed, 1 skipped`. PLUS a FULL-suite run (deselect #212,
   `-k` not sphere apply-bit-id): **0 EmissionSpectrum reds** (definitive —
   no runtime-computed-chi site escaped the AST scanner).
4. **pyright:** `2282 errors, 19 warnings` — **EXACTLY baseline 2282, 0
   net-new**. The `Mixture.is_fissile` skeleton diagnostic cleared (property
   landed). **0 new `# type: ignore`** in all 3 production files. The 4
   focused errors in mixture.py/isotope.py are ALL PRE-EXISTING (compute_macro_xs
   `int|ndarray` SigP/Sig2, Isotope.is_fissile np.bool_ return at :71) — the
   global count is unchanged so none are net-new. (The pyright RATCHET test
   fails because counts DECREASED codebase-wide — `data:59→8`, `sn:263→173`,
   etc.; IDENTICAL at clean HEAD via git stash, a PRE-EXISTING un-re-baselined
   condition from prior #226 cleanup, NOT S10a — NO module INCREASED.)

**Sphinx:** `build succeeded, 1 warning` — the 1 warning is a PRE-EXISTING
Python SyntaxWarning (unescaped LaTeX in `test_projection_operators.py`/
`test_fission_operator.py` docstrings, captured by the Nexus graph-builder;
files I did NOT touch; IDENTICAL count at clean HEAD via git stash). My new
section/label/todo add ZERO new warnings (no undefined label / unknown ref).
Label `emission-spectrum-simplex-law` rendered in the built HTML.

## Self-improvement / notes

- NO new vv failure mode (gates are vv #11 contract-validation both-legs +
  Mode-8, both in the skill). NO ERR (S10a is ADDITIVE; the dead non-fissile
  χ was inert-not-wrong; the SigF/SigP seam is a design subtlety not a caught
  bug). ERR-063 stays reserved.
- **LESSON (reinforces L20):** a "guard-at-the-data-source" change has a blast
  radius = EVERY direct constructor of the guarded type, NOT just the factory
  path. The factory-path audit (`get_mixture` callers) is necessary but NOT
  sufficient — DIRECT `Mixture(...)`/`make_mixture(...)` callers with inline
  placeholder fields are invisible to it. A deterministic AST scanner (find
  every call to the constructor, statically classify the guard-keying field +
  the guarded field) found all 10 main-tree sites in one pass where running the
  suite suite-by-suite leaked them one error at a time. **A `bool()`-wrapped
  predicate (`Mixture.is_fissile`) is pyright-clean; the un-wrapped twin
  (`Isotope.is_fissile`, `np.bool_`) is a pre-existing red — mirror the wrap
  when adding the sibling.**
- **LESSON (new seam):** a 2-field fissility signal (`SigF` the cross-section
  vs `SigP=νΣf` the production rate) can DISAGREE in synthetic fixtures that
  set production directly. A guard keyed on ONE of them (`SigF`) can demand
  `χ=0` on a material that genuinely emits via the OTHER (`SigP>0`). The
  behavior-preserving fix is to make the FIXTURE consistent (set the keyed
  field), NOT to zero the consumed χ. Flag the keying choice to the user — it
  is a real design question deferred past S10a.

Extends [[issue-257-s5-functional-category-closeout]] (same campaign, the
typed-field-algebra family) + [[issue-257-s10a-emission-spectrum-verification]]
(the TA gate-spec this implements).

## ⭐ S10a REFINEMENT (2026-06-21, review outcome — elegance+qa) — guard re-keyed on PRODUCTION

Applied the review outcome on the SAME branch `feature/field-typed-operator-algebra`
(NOT committed). The S10a guard was re-keyed from FISSIONABILITY (`SigF>0`)
to PRODUCTION (`νΣf>0`), which BOTH retires the SigF/SigP hack above AND is
physically correct: χ is consumed ONLY as a fission SOURCE `χ·νΣf·φ`, so a
valid emission spectrum is required exactly where production is nonzero.

**Production edits (3 files):**
1. `emission_spectrum.py` — NEW module-level helper `enforce_emission_spectrum(chi, *, is_producing)` that single-sources the coerce-and-branch law (Pattern 2 — the two `__post_init__` bodies were duplicated): `chi=EmissionSpectrum(chi); chi.assert_normalized() if is_producing else chi.assert_null(); return chi`. Also sharpened the `is_emitting` docstring to make it a SPECTRUM-level query ("non-null, any χ_g>0") DISTINCT from a material's production predicate. `is_emitting` KEPT (rounds out the value-object's natural API).
2. `isotope.py` — NEW `is_producing` property `bool(np.any(self.nubar*self.sigF>0))` (broadcast: `sigF` `(n_sig0,NG)`, `nubar` `(NG,)` → `(n_sig0,NG)`). `__post_init__` routes through the helper. **`is_fissile` (`np.any(sigF>0)`) KEPT UNCHANGED** — it is the distinct "can fission?" predicate consumed by `compute_macro_xs` `fissile_indices`; sharpened its docstring to contrast with `is_producing`.
3. `mixture.py` — the S10a-added `is_fissile` property (Mixture had none before; grep-confirmed ZERO production consumers, only the gate test) REPLACED by `is_producing` `bool(np.any(self.SigP>0))` (SigP=νΣf). `__post_init__` routes through the helper.

**Hack RETIRED:** `test_trajectory_resolvent_billiard.py::_mixture_from_xs` — the `SigF=nu_sf_arr.copy()` stand-in REVERTED to `SigF=np.zeros(ng)`. The material is PRODUCING (`SigP=νΣf>0`), so `is_producing` is True and its default-simplex χ is correctly required — no SigF stand-in needed (the billiard path reads SigP/χ, never SigF). The whole SigF/SigP "design seam" from the original closeout DISSOLVES: production-keying classifies the production-bearing fixture correctly with no special handling.

**Gate test re-keyed** (`test_chi_invariant_enforcement.py`): `_mixture(*, sig_p, chi)` (parametrize by SigP); `_isotope(*, sig_f_row, nubar, chi)` (producing legs set nubar>0 where sig_f_row>0). Methods renamed producing/non-producing; assert `is_producing`. NEW `test_non_producing_nonzero_raises` ISOLATES the predicate gap: sigF>0 (fissile) but nubar=0 (non-producing) → χ must be null → proves `is_producing` ≠ `is_fissile`. Real-GENDF leg asserts `is_producing` (U_235 True, O_016/H_001 False — VERIFIED on real `.h5`: U_235 max(νΣf)=3.17e4, O_016/H_001 ≡0). BOTH legs/branches kept (vv #11).

**Promoted ERR-063 test:** `derivations/diagnostics/diag_err063_probe_e_intrinsic_property.py` → NEW `tests/derivations/test_peierls_fission_source_indexing.py` (`@l1 @catches("ERR-063") @verifies("peierls-mg-operator")`, Mode-8-safe `np.testing.assert_array_less`). Claim: a non-fissile region's χ must not affect k_eff (Hébert 2009 Eq.3.57/3.58). MUTATION-VERIFIED teeth on the FINAL file: reverting `geometry.py:6426` to `chi_n[i,ge]` reddens BOTH negative legs (3.59% k_eff move) while the positive non-degeneracy leg stays green; restored via Edit (NEVER `git checkout` an uncommitted tracked file — elegance lesson). ⚠ PYRIGHT TRAP: the diagnostic lived in `derivations/diagnostics/` (UNSCANNED) and used `solve_peierls_mg(geom, **case_dict, **_QUAD)` `**dict` splat → the sibling peierls MG tests carry 66 `reportArgumentType` from this idiom, but promoting into `tests/derivations/` (SCANNED) added **+41 net-new pyright**. RESOLVED to 0 net-new (NO `# type:ignore`) by collapsing both call sites into ONE typed helper `_solve_two_region(geometry, *, chi)` with FULLY-EXPLICIT kwargs (no splat) + a `pytest.fail` None-narrow on `k_eff: float|None` (Mode-8-safe).

**Verification (all `-O`, route around #250/#232/#212):**
- `test_emission_spectrum.py`+`test_chi_invariant_enforcement.py`: **28 passed**.
- DD regression byte-id: **13 passed** (13 within-tol pre-existing DriftWarnings, NO snapshot moved).
- promoted test + rank_n_class_b_mr_mg + billiard: **42 passed, 1 skipped, 2 xfailed** (467s; billiard slow).
- TA full-regression command: **107 passed, 2 xfailed** (pre-existing operator-subset).
- ALL 11 edited-site suites (data/moc/cross_method/scattering/axis/d3/sweep): **356 passed, 1 skipped, 0 failures** under production-keying.
- promoted test teeth: reverted-`chi[i]` → **2 failed (negative legs RED) + 1 passed (positive green)**, restored → 3 passed.
- **pyright: 2353 errors EXACTLY = true baseline** (measured via git-stash of tracked edits + moving untracked S10a files aside; the closeout's earlier 2282 was STALE branch-drift — the real pre-S10a count is 2353). **0 net-new, 0 new `# type:ignore`.**
- Sphinx `-W`: **build succeeded** (`enforce_emission_spectrum` `:func:` resolves; `peierls-mg-operator` verifies-edge points at the now-CORRECTED source-indexed equation block; the LD `verifies(...) — skipping` notices are pre-existing graph-builder INFO, NOT my refs).

**Sphinx stub UPDATED** (`cross_section_data.rst`, STILL a stub — rich narrative owed to archivist): the `emission-spectrum-simplex` label math + prose re-keyed fissile→producing; `enforce_emission_spectrum` `:func:` added; TODO point (5) rewritten from "the SigF/SigP subtlety" to "why the law keys on PRODUCTION not fissionability + the two-predicate divergence on synthetic SigP>0/SigF=0 fixtures".

**NOT touched** (per brief): peierls `geometry.py`/`slab.py` source-index fix (untouched except the temp teeth-probe, fully restored byte-identical), `peierls_nystrom.rst` theory docs (already corrected), the xs_library non-fissile χ zeroing + the 10 AST-found precursor fixtures (now inert under BOTH the peierls fix AND production-keying), NO χ numerical value changed, NO commit.

**LESSON (the SigF/SigP seam RESOLUTION, supersedes the original closeout's "flag to user"):** when a guard-at-source keys on a predicate, key it on the predicate that matches WHY the guarded value exists — χ exists to be a fission SOURCE (`χ·νΣf·φ`), so PRODUCTION (νΣf>0) is the right key, NOT fissionability (Σf>0). The two coincide for real `compute_macro_xs` data (νΣf>0 ⟺ Σf>0) but diverge for synthetic production-bearing fixtures (SigP>0, SigF=0); production-keying classifies them correctly with NO stand-in, dissolving what looked like an irreducible "design seam". The fix retires a hack instead of documenting one.
**LESSON (pyright scope trap on diagnostic promotion):** promoting a diagnostic from an UNSCANNED dir (`derivations/diagnostics/`) into a SCANNED test dir (`tests/`) imports its pyright debt — a `**dict`-splat into a typed solver signature is invisible noise in the diagnostic but +N net-new once scanned. Collapse the call to a single typed helper with explicit kwargs (no splat) to clear it without `# type:ignore`; narrow a `float|None` return with `pytest.fail` (Mode-8-safe), not bare `assert`.
