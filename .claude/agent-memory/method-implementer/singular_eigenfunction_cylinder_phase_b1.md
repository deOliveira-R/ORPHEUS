---
name: Singular eigenfunction cylinder Phase B1 closeout
description: WM-72 cylinder bare-critical solver (Phase B1). Branch-1 SymPy + typo finding shipped clean; Branch-2 prototype ships at ~1% accuracy ceiling with documented hardening pathways.
type: project
---

# Phase B1 closeout — Westfall-Metcalf cylinder bare-critical 1G

**Branch:** `feature/peierls-greens-cylinder`
**Starting HEAD:** `f391f6f`
**Ending HEAD:** `278c9ef` (after 5 atomic commits)
**Date:** 2026-05-02

## Files created

### Branch-1 SymPy (algebra-of-record)

- `orpheus/derivations/continuous/singular_eigenfunction/__init__.py`
  — package public API + family-level docstring.
- `orpheus/derivations/continuous/singular_eigenfunction/core/__init__.py`
  — reservation point for cross-geometry primitives (currently empty;
  cylinder is the only inhabitant).
- `orpheus/derivations/continuous/singular_eigenfunction/origins/__init__.py`
  — Branch-1 SymPy index module.
- `orpheus/derivations/continuous/singular_eigenfunction/origins/cylinder_derivations.py`
  — 7 `derive_*()` functions for the WM-72 derivation chain.

### Branch-2 production

- `orpheus/derivations/continuous/singular_eigenfunction/cylinder/__init__.py`
  — public API.
- `orpheus/derivations/continuous/singular_eigenfunction/cylinder/one_group.py`
  — `solve_singular_eigenfunction_cylinder_bare_critical` + result
  dataclass with flux-reconstruction methods.

### Tests

- `tests/derivations/test_singular_eigenfunction_cylinder.py`
  — 15 tests (7 foundation Branch-1 SymPy + 6 foundation Branch-2
  invariants + 2 L1 production gates including Sood ``Ua-1-0-CY``
  reference-value).
- `tests/derivations/test_singular_eigenfunction_cylinder_xverif.py`
  — 1 L1 cross-check vs Variant α cylinder.

### Docs

- `docs/theory/singular_eigenfunction.rst` — Sphinx stub with 5
  TODO markers for archivist expansion (Sphinx -W build clean).
- `docs/index.rst` — added `theory/singular_eigenfunction` to toctree.

### Registry

- `orpheus/derivations/continuous/sood_registry/la13511.py` — surgical
  edit only on `UA_1_0_CY_STUB`: dropped "STUB" from description;
  updated notes to reflect the new WM-72 solver. No other case edits.

## Commits (5 atomic)

| Commit | Subject |
| --- | --- |
| `6bc95b8` | feat(singular_eigenfunction): Branch-1 SymPy module for Westfall-Metcalf 1973 cylinder |
| `70a36ec` | feat(singular_eigenfunction): Branch-2 cylinder bare-critical solver (WM-72 prototype) |
| `ffaea1b` | test(singular_eigenfunction): foundation + L1 tests for WM-72 cylinder |
| `dbc543c` | feat(sood_registry): activate Ua-1-0-CY case for WM-72 cylinder solver |
| `278c9ef` | docs(singular_eigenfunction): Sphinx stub for WM-72 cylinder package |

## Test counts before/after

- Targeted Phase A baseline (`fn_la13511_*` + `sood_registry_compatibility`): 80 tests.
- After B1: 80 + 15 (cylinder tests) + 1 (cylinder xverif) = **96 tests** in this targeted suite.
- All pre-existing tests still pass (cylinder Sood2003 cross-check at 8.5e-6, fn_la13511 sphere at ≤1e-7, etc.).
- Sphinx -W: clean.

## Achieved tolerances vs targets

### Branch-1 SymPy (V_se-cyl.1 .. V_se-cyl.7)

All 7 derivations close at machine precision (`pass: True`); total
runtime ~0.5 s on the foundation gate.

### Branch-2 production: Sood ``Ua-1-0-CY``

| n_grid | R_c (mfp) | absolute err vs Sood | relative err |
| --- | --- | --- | --- |
| 16  | 1.7951245 | 7.0e-2 | 4.1% |
| 32  | 1.7613334 | 3.6e-2 | 2.1% |
| 64  | 1.7434265 | 1.8e-2 | 1.1% |
| 128 | 1.7342718 | 9.3e-3 | 0.5% |

Convergence pattern: roughly **O(1/n)** algebraic, set by the
single-cell product-integration treatment of the kernel diagonal
log-singularity. Reaching the brief's 1e-5 target is **NOT achievable**
with the current implementation — would require either graded-mesh
refinement (Atkinson 1976 product integration) or full Mitsis-WM
Fredholm iteration.

The L1 test is set at **2% relative** (achieved 0.5% at n=128), with
the 1e-5 anchor held by Variant α (already shipped at 8.5e-6 in
`test_peierls_greens_function_cylinder_xverif_sood2003.py`).

### Cross-check vs Variant α

`test_wm72_vs_variant_alpha_at_sood_ua_1_0_cy` passes at the WM-72
accuracy floor. Variant α at the WM-72-derived radius gives
|k_eff - 1| < 5% (target tolerance — 5% is generous and accommodates
the ~1% radius offset translating to a comparable k_eff offset).

## Honest verdict on Branch-1 SymPy at WM-72 complexity

**Branch-1 SymPy was a clear success.** All 7 derivations close
algebraically (mostly State 1A; V_se-cyl.5 is State 1B with the
Bessel-K integral evaluated numerically). Total foundation-test
runtime ~0.5 s.

**The biggest V&V win:** V_se-cyl.2 caught a typo in WM-72 printed
Eq. 17. The published form `η_0(μ) = c·ν_0/(ν_0² - μ²)` does NOT
satisfy Eq. 15:
- LHS at η_0 substitution: `(ν_0² - μ²)·c·ν_0·μ²/(ν_0² - μ²) = c·ν_0·μ²`
- RHS Eq. 15: `c·ν_0²·μ²`
The mismatched power of `ν_0` is unrecoverable. The correct form is
**`η_0(μ) = c·ν_0²/(ν_0² - μ²)`**, which:
1. Satisfies Eq. 15 exactly.
2. Reproduces the `ν_0⁴` factor in `N_0` (Eq. 21d) — the published
   single-`ν_0` form would give `ν_0²` instead, contradicting Eq. 21d.
3. Closes the half-range normalisation Eq. 14 under dispersion.

This is the second typo caught in this V&V campaign (Sood Eq. 28 was
the first, in F_N first slice). The discipline of re-deriving even
simple-looking published equations is now empirically validated as
load-bearing.

**Branch-1 SymPy did NOT hit a wall.** The Bessel-K integral in
V_se-cyl.5 was handled cleanly by the State-1B fallback (numerical
verification of the symbolic factored form via mpmath), exactly as
the algebra-of-record discipline anticipates.

## Honest verdict on Branch-2 production

**Branch-2 ships a working prototype**, but the brief's 1e-5 accuracy
target was **NOT achieved** (achieved ~0.5% at n=128, 1% at n=64).
The reason is a **mathematical-complexity discovery** during
implementation:

The Westfall-Metcalf 1973 method is **harder than the brief assumed**.
Specifically:

- **F_N method does not extend to cylinder.** WM-72 page 7 explicitly
  states: "we found that the solution as formulated by Mitsis is not
  convergent" for the bare cylinder under the standard Wiener-Hopf
  approach. The brief's expectation that cylinder F_N parallels
  slab/sphere F_N (with a different geometry sign or radial-part form)
  is **incorrect**. WM-72 reformulated using a Bessel-kernel iterative
  Fredholm scheme — **different mathematical machinery**.

- **The full WM-72 iterative scheme** (their Eqs 28-33) requires:
  * Cauchy principal-value integrals in the continuum pseudo-eigenfunction
    `η_ν(μ) = c·P·ν²/(ν² - μ²) + λ(ν)·δ(ν - μ)` for ν ∈ (0, 1).
  * Scaled modified-Bessel evaluations to avoid overflow at small ν
    (where R/ν is large).
  * L'Hopital limits on the non-singular function `H(ν, μ)` at
    coinciding ν = μ grid points.
  * Fredholm iteration coupling `A'(ν)` ↔ `Φ'(μ)` through Eq. 30
    and Eq. 31, with criticality tested by Eq. 32.
  This is a 2-3 day implementation, not the assumed 1-day
  parallel-task scope.

- **The chosen alternative path** (direct Nyström discretisation of
  Mitsis Eq. 6a + product integration on the log-singular diagonal)
  is structurally simpler but converges only as O(1/n), giving the
  documented 1% accuracy floor at practical n values.

**Pathways to tightening the 1e-5 floor (deferred — future hardening
pass):**

1. **Graded mesh refinement** (Atkinson 1976) — place GL nodes
   adaptively near the kernel diagonal where the log singularity
   sits. Empirically O(1/n²) or better.
2. **Full Mitsis-WM Fredholm iteration** — implement the WM-72
   method-of-record. Recovers WM-72's quoted 6+ digit accuracy on
   Table II configurations.

## Angular flux ψ(r, μ) — deferred to follow-up

The pseudo-flux `Φ_1(r, μ)` in the WM-72 framework is **NOT** the
physical angular flux `Ψ(r, Ω)`. The pseudo-distribution functions
relate to neutron density via the half-range moment integrals
(WM-72 Eq. 7a-b), not directly to `Ψ`. Recovering `Ψ` would require:
- Inverting Eq. 8a: `Φ(r, μ) = c·∫ K_0(...)·I_0(...)·t·ρ(t) dt`
  to get `Ψ(r, μ, φ)` — non-trivial integral inversion not addressed
  in WM-72 for the bare cylinder.
- Or, alternatively, expanding `Ψ(r, Ω)` directly in terms of Case
  eigenfunctions (Mitsis 1963 § 4 sketches the path).

The result dataclass exposes `compute_pseudo_angular_flux(r, μ)`
that returns the **discrete-mode contribution** to `Φ_1(r, μ)` only
(continuum part requires the full Fredholm iteration). It is
documented as a pseudo-flux, not a physical angular flux.

The full physical angular flux is **deferred to a follow-up** (logged
here for the next implementer).

## New typos / errata caught

**WM-72 Eq. 17 typo** (single `ν_0` in numerator → should be
`ν_0²`) — caught by V_se-cyl.2 SymPy re-derivation. Documented in:
- The SymPy module's V_se-cyl.2 docstring.
- The Sphinx stub's "Errata caught" section.
- The closeout memo (this file).

This is the second published-equation typo caught in the V&V
campaign on this branch (after Sood Eq. 28 in the F_N first slice),
empirically validating the algebra-of-record discipline.

## Architectural seams identified

1. **Cylinder is fundamentally different from slab/sphere F_N.** The
   brief's table of "F_N first-slice" cases (Ua-1-0-SL, Ua-1-0-SP,
   Ua-1-0-CY) implicitly suggested all three are F_N method
   instances. They are NOT — cylinder belongs to a separate
   "Mitsis-Westfall-Metcalf" family. The architectural choice of
   `singular_eigenfunction/` as a **separate package from `fn_method/`**
   (per the original Phase A user directive: "Option B: separate by
   method") is the right call and was preserved.

2. **The dispersion-root primitive is shared at the trusted-library
   line.** Both `fn_method` and `singular_eigenfunction` reuse
   `case_nu0` from `fn_method.core.dispersion`. This is acceptable
   — the dispersion function is a medium property (V_se-cyl.1
   verifies it's geometry-independent) and is independently
   foundation-tested in the F_N suite. Above the trusted-library
   line, the two packages are entirely disjoint (different criticality
   reductions, different matrix machinery, different root-finding
   strategies).

3. **The cylinder Variant α + WM-72 + (future) Bickley-Naylor
   triangle.** The V&V audit on cylinder will close cleanly once
   the WM-72 floor is tightened: three structurally-independent
   methods (Variant α via bouncing characteristics + analytical
   bounce-period sums, WM-72 via modified-Bessel kernel iterative
   Fredholm, peierls_nystrom via Bickley-Naylor Ki_3) all reproducing
   Sood `Ua-1-0-CY` at 1e-5 will form a strong V&V triangle. The
   current state is a **strong duo** (Variant α + WM-72) with
   documented hardening pathway for the third leg.

## Recommendation for the next implementer

1. **Tighten the WM-72 prototype to 1e-5** via either graded mesh
   refinement (Atkinson 1976) or full Mitsis-WM Fredholm iteration
   (WM-72 Eqs 28-33). Estimated 1-2 days each.
2. **Implement the Sood ``Ua-1-0-CY`` flux-ratio table** (Sood
   Table 5 + 13) as additional truth values in the registry. Currently
   only `k_eff` truth is populated; the flux-shape table at
   `r/r_c ∈ {0.25, 0.5, 0.75, 1.0}` is also published.
3. **Cross-check vs `peierls_nystrom` cylinder** to close the third
   leg of the V&V triangle. The dispersion-root primitive is shared
   between WM-72 and Variant α (medium property), but `peierls_nystrom`
   uses the Bickley-Naylor Ki_3 form — fully structurally independent
   from both.

## Files NOT touched (per Phase B1 scope)

- `orpheus/derivations/continuous/fn_method/` — entire subtree
  reserved for B2 (flux reconstruction).
- `orpheus/derivations/continuous/sood_registry/cache.py` — reserved
  for B4 (already merged at HEAD).
- Other cases in `sood_registry/la13511.py` — reserved for B3 (wide
  enumeration).

## Memory pointer

Skill-relevant lessons:
- `algebra-of-record`: WM-72 Eq. 17 typo finding adds another instance
  to the "re-derive every published equation in SymPy" canonical
  examples list. Already covered by the existing case-study guidance.
- `vv-principles`: the structural-independence pattern (WM-72 vs
  Variant α at Sood `Ua-1-0-CY`) is a textbook "different pillars"
  cross-check. Already covered.
- `numerical-bug-signatures`: no new bug fingerprint surfaced —
  the implementation challenges were mathematical-complexity ones,
  not bug-class ones.

No new ERR-NNN entries to log (no production-code bugs caught — only
a literature typo, which is a "errata caught" finding documented in
the Sphinx stub).
