---
name: issue-257-s10a-emission-spectrum-verification
description: S10a gate spec — EmissionSpectrum(ndarray) newtype + Mixture/Isotope __post_init__ simplex/null χ guard; CRITICAL precursor finding (377-symbol blast radius from non-fissile xs_library chi); byte-id anchors confirmed.
metadata:
  type: project
---

# #257 S10a — EmissionSpectrum newtype + container χ-invariant guard (PRE-IMPL gate spec)

Branch `feature/field-typed-operator-algebra`. SUT = a self-validating
value-object for the fission emission spectrum χ, enforced at the data
source (Mixture + Isotope `__post_init__`). S10a = TYPE + GUARD only;
χ-value change is the deferred S10b (multi-fissile production-weighting).
NO production code written by this agent — spec + skeletons only.

## ⭐⭐ CRITICAL PRECURSOR FINDING — the brief's L20 audit is INCOMPLETE

The brief states "21 Mixture construction sites; 19 comply; exactly 2
need a precursor edit (`test_sweep_cache.py:287`/`:598`, confirmed at
those lines)". That count covers only DIRECT `Mixture(...)`
constructions in test bodies. It MISSES the dominant factory path:
**`orpheus.derivations.common.xs_library.make_mixture` / `get_mixture`**.

PROBED at HEAD c6e21c0:
- `make_mixture` does `chi=chi.copy()` from the region dict.
- The non-fissile regions B, C, D hardcode NON-ZERO chi:
  `_B_1G/_C_1G/_D_1G` carry `chi=np.array([1.0])`; the 2G variants
  `chi=[1.00, 0.00]`; the 4G variants `chi=[0.60,0.35,0.05,0.00]`.
- `get_mixture("B"/"C"/"D", ng)` → a non-fissile `Mixture`
  (`SigF≡0`) whose NEW `__post_init__` would call `assert_null()` and
  RAISE on the non-zero chi.
- Blast radius (Nexus `impact` upstream on `get_mixture`):
  **377 total affected symbols**; **77 non-fissile
  `get_mixture("B"/"C"/"D", ...)` call sites in `tests/`** — SN, CP,
  MoC, MC verification suites + the SN regression-snapshot generator
  (`tests/sn/regression/_generate_snapshots.py` uses `get_mixture("B",
  ng)` pervasively). ALL would red en masse the moment the non-fissile
  branch is enabled.

⇒ This is the canonical the L20 retirement-scope hazard (factory-path callers dwarf direct-construction callers) retirement/scope-blowup hazard
manifesting at design time. **The precursor is NOT 2 edits; it is:**
1. The 2 direct-`Mixture` fixtures (`test_sweep_cache.py:287/:598` →
   `chi=np.zeros(1)`), AND
2. **Zero the non-fissile chi in `xs_library.py` regions B/C/D** (the
   `_B_*`/`_C_*`/`_D_*` dicts: `chi=np.zeros(ng)`). This is
   BEHAVIOR-NEUTRAL: for non-fissile mixtures `SigP=νΣf≡0`, and
   `FissionOperator.apply` = `χ·(νΣf·φ).sum(...)` → the χ value is
   multiplied by zero in EVERY fission/keff path. Zeroing χ cannot move
   any converged flux or k_eff (proven by the byte-id anchor below).

The method-implementer + main agent MUST land BOTH precursor parts
BEFORE enabling the `__post_init__` non-fissile branch. Recommend the
main agent dispatch `explorer` for a full production-caller audit of
`make_mixture`/`get_mixture` (Sood registry, peierls_nystrom cases,
diffusion/cp solvers also call it) — there may be NON-test production
callers (e.g. `derive_sn_heterogeneous_continuous`, peierls cases) that
also build non-fissile mixtures with the library chi. The 77 figure is
tests-only; the 377 is the full graph.

## SUT shape (recap)

- `EmissionSpectrum(np.ndarray)` in NEW leaf
  `orpheus/data/emission_spectrum.py` (cycle-safe). Methods operate on
  the array's OWN values: `assert_normalized()` (Σχ≈1 AND χ≥0),
  `is_emitting` property (`bool(np.any(self>0))` — Python bool by spec),
  `assert_null()` (χ≡0 strict). Needs `__array_finalize__` (flagged for
  implementer).
- `Mixture` (`macro_xs/mixture.py`, @dataclass, NO `__post_init__`
  today) + `Isotope` (`micro_xs/isotope.py`, @dataclass, NO
  `__post_init__`, ALREADY has `is_fissile = np.any(self.sigF>0)`) gain
  `__post_init__`: coerce `self.chi`→EmissionSpectrum, then conditional
  `if is_fissile: assert_normalized() else: assert_null()`. Mixture
  needs a NEW `is_fissile` property mirroring Isotope's (`np.any(SigF>0)`
  — returns `np.bool_`, NOT Python bool; see teeth note below).
- χ flow: `Isotope.chi (ng,)` → `Mixture.chi (ng,)` → per-cell
  `MaterialXSField._chi_cell`/`FissionOperator.chi` → `RankOneOperator`.
  Invariant asserted ONLY at the source — NEVER per-cell (non-fissile
  cells legitimately hold χ=0).
- FP tol (probed): real `load_isotope('U_235',294)` → `chi.sum()==
  1.0000000000000002` (residual 2.22e-16). Recommend
  `np.isclose(self.sum(), 1.0, atol=1e-12, rtol=0)` (~4 orders headroom).

## Gates (test skeletons LANDED, pre-impl)

Two files under `tests/data/`, both `@pytest.mark.foundation` (software
invariant — probability-simplex law has NO theory `:label:`; foundation
carries NO `verifies(...)` per [[feedback_vv_tagging]]):

### `tests/data/test_emission_spectrum.py` — gate 1 + gate 3a
- Whole module `pytest.importorskip("orpheus.data.emission_spectrum")`
  → COLLECTS-SKIPPED pre-impl, green post-impl, ZERO edits needed.
- `TestNdarraySubclassBehaviour` (gate 3a, zero-ripple at type level):
  ndarray-subclass, `np.asarray(ES(a))==a` roundtrip, `.sum()`,
  `chi[None,:]` broadcast, `.copy()`, einsum-feed value — pins the
  `__array_finalize__` contract that lets existing call sites work
  unchanged.
- `TestAssertNormalized` (vv #11 both legs): simplex `[.6,.35,.05]`
  PASSES; `[.6,.35,.10]` (Σ=1.05) RAISES (sum clause); `[1.2,-0.2]`
  (Σ=1 but negative) RAISES (negativity clause — SEPARATE from sum,
  the clause-isolation leg); tol edges `[0.5, 0.5+2e-16]` PASSES,
  `[0.5, 0.5+1e-6]` RAISES (pins both ends of the 1e-12 band).
- `TestIsEmitting` / `TestAssertNull` both legs.

### `tests/data/test_chi_invariant_enforcement.py` — gate 2 + gate 5
- HAND-BUILT minimal fixtures `_mixture`/`_isotope` (L11 — NEVER via
  `make_mixture`/`load_isotope` for the synthetic legs; those builders
  ARE the guarded system).
- `TestMixturePostInit` + `TestIsotopePostInit` (gate 2, vv #11, BOTH
  branches BOTH legs): fissile+simplex→ok; fissile+non-simplex→raises;
  non-fissile+zeros→ok; non-fissile+nonzero→raises; chi-coerced-to-ES.
  The `non_fissile_nonzero_raises` leg's docstring names the 77-site
  precursor hazard.
- `TestRealGendfConstructs` (gate 5, the production-path no-red proof):
  uses the REAL `load_isotope` (deliberately — the flip side of L11) —
  `U_235` fissile constructs w/o raising + residual inside band;
  `O_016`/`H_001` non-fissile (χ≡0 on real data) construct w/o raising.
  `skipif` guarded on `.h5` presence.

### Mode-8 (gate 4) — TEETH-dry-run-PROVEN
- ORPHEUS runs `python -O` → bare `assert` STRIPPED. ALL substantive
  legs route through `pytest.raises` / `np.testing.*` / a module-level
  `_require(cond,msg)` helper (`pytest.fail` — a CALL, fires under -O).
  ⚠ CAUGHT IN DRY-RUN: my first draft used `assert bool(m.is_fissile)`
  / `assert isinstance(...)` positive legs — STRIPPED under -O (the
  `PytestConfigWarning` flagged it). Converted ALL to `_require(...)`.
  Also caught: `is_fissile` returns `np.bool_` not Python `bool` →
  `assert ... is True` would FALSELY RED post-impl → use `bool(...)`
  truthiness. (`is_emitting` IS Python bool by spec — `is True` kept,
  + a `type(...) is bool` leg to pin the `bool(...)` coercion.)
- PRE-IMPL `-O` baseline (PROVEN, honest): container file = **6 failed
  / 5 passed / 2 skipped**. The 6 fails are the RIGHT pre-impl signals:
  2 Mixture positive legs = `AttributeError: no is_fissile` (unblocks
  when the property lands); 4 negative legs = `DID NOT RAISE` (unblock
  when `__post_init__` lands). The 5 passes = Isotope positive legs
  (is_fissile already exists) + real-GENDF legs (real chi already
  valid). 2 skips = ES-coercion legs (importorskip). emission_spectrum
  file = 1 skipped (whole-module importorskip). This is teeth-proven:
  the negative legs RED for the absence of the guard and flip green
  ONLY when the guard lands.

## Gate 3 — byte-identity anchors (CONFIRMED GREEN at HEAD c6e21c0)

The headline claim is "χ VALUES unchanged, type wraps them". Two anchor
tiers, both PROBED passing at HEAD:
1. **Strictest (snapshot byte-id):**
   `tests/sn/regression/test_dd_regression.py` — 13 passed (45.8s; the
   13 DriftWarnings are WITHIN-tol, pre-existing, NOT a S10a effect).
   Its `_generate_snapshots.py` builds `get_mixture("B", ng)` non-fissile
   mixtures pervasively → AFTER the precursor zeroes xs_library chi,
   these snapshots MUST stay byte-identical. THE proof that non-fissile
   chi is inert in the converged flux. ⚠ do NOT regenerate snapshots
   for S10a (zeroing inert chi changes no value — if a snapshot moves,
   STOP: it means chi was NOT inert and the precursor is unsafe).
2. **FissionOperator / keff path (named in brief):**
   `test_fission_operator.py::TestBitIdenticalExtraction`
   (`get_mixture("A","2g")` fissile simplex + `get_mixture("B","2g")`
   non-fissile; `nulp=4`) + `test_fission_kernel_crosscheck.py` +
   `test_kinf_homogeneous.py` (L1 `A⁻¹F` analytical, struct-indep) +
   `test_keff_2d.py::...test_si_krylov_heterogeneous_2g_nonflat_flux` —
   ALL passed at HEAD (50 passed / 2 xfailed pre-existing for the
   operator subset; keff_2d 1 passed).
- Data-layer green-anchors (cheap, must stay green): `test_mixture.py`,
  `test_cross_section_data.py` (both passed).

## Regression command (route around baseline reds #250/#232/#212)

Run AFTER precursor + production land, under `-O` (canonical):

```
.venv/bin/python -O -m pytest \
  tests/data/test_emission_spectrum.py \
  tests/data/test_chi_invariant_enforcement.py \
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
  w/ #232) if any surface in the regression run.
- DESELECT `tests/sn/eigenvalue/test_keff_slab.py::...heterogeneous_absolute_keff`
  (#212 `continuous_get` hang) if you widen to the keff suite —
  `--deselect tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff`.
- These reds are PRE-EXISTING on the branch (NOT S10a) — confirm at
  clean HEAD before crediting any S10a red to S10a.

## S10b SEAM (where production-weighted χ_mix verification attaches)

S10a leaves `compute_macro_xs` (mixture.py:191-194) using the FIRST
fissile isotope's χ verbatim (`chi = isotopes[fissile_indices[0]].chi.copy()`)
— a documented "simplification". S10b replaces this with the
production-rate-weighted mixture spectrum
`χ_mix,g = Σ_i N_i νΣf_i χ_i,g / Σ_i N_i (νΣf_i·1)` (multi-fissile
weighting). The S10a guard is the SCAFFOLD that S10b plugs into:
- S10b's new `χ_mix` is a value-CHANGE → it will be a NON-bit-identical
  regression (the [[feedback_principled_over_bit_identical]] gate, 3
  criteria; struct-indep ref = hand-weighted χ_mix vs a 2-fissile-isotope
  mixture, NOT another solver).
- S10a's `EmissionSpectrum.assert_normalized()` is REUSED VERBATIM as
  S10b's correctness floor: a correctly production-weighted convex
  combination of simplices is itself a simplex → the S10a guard fires
  automatically on the S10b output (a weighting bug that breaks the
  simplex is caught for free at construction). S10b's NEW gate is the
  VALUE (the weights), not the simplex-ness.
- S10b attaches at `compute_macro_xs` (the producer) — add a
  `test_chi_mix_production_weighting.py` with a 2-fissile-isotope
  Mixture whose hand-computed χ_mix differs from either input χ;
  S10a's container guard is the no-extra-work invariant leg.

## Self-improvement
- NO new vv failure mode (gates are vv #11 contract-validation + Mode-8,
  both already in the skill table). NO ERR (next free ERR-063 — no real
  bug, S10a is additive).
- The PRECURSOR finding is a the L20 retirement-scope hazard (factory-path callers dwarf direct-construction callers) reinforcement (factory-path
  callers dwarf direct-construction callers in a guard-at-source change)
  — captured here, not a skill edit (skill stays uncommitted per L28/
  cluster policy; main authorizes skill edits).

Extends [[issue-257-s5-functional-category-verification]] (the same
campaign, the typed-field-algebra family) + [[feedback_test_intrinsic_properties]]
(intrinsic-law gate per math-bearing type) + [[feedback_vv_tagging]]
(foundation NO verifies) + [[feedback_regression_tolerance_design]]
(byte-id vs principled, the S10b seam) + the L20 retirement-scope hazard (factory-path callers dwarf direct-construction callers) (the precursor
scope-blowup).
