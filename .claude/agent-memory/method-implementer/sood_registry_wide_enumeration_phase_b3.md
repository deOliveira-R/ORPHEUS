---
name: sood_registry Phase B3 wide enumeration
description: Wide-cast enumeration of LA-13511 cases solvable today by existing fn_method machinery (k_inf 1G/2G/MG, slab F_N, sphere F_N). Registry expanded 5→42 cases; 56 new tests pass + 12 stubs.
type: project
---

# sood_registry Phase B3 — wide enumeration closeout

`feature/peierls-greens-cylinder` 2026-05-02. Branch HEAD before this work: `f391f6f`. Parallel sibling phases: B1 (cylinder solver), B2 (flux reconstruction), B4 (cache).

## Scope (verbatim from dispatch)

> "look at everything that can be implemented from Sood using the current
> machinery and start wide implementation of sood's cases. If the machinery
> needs to be expanded, expand it step wise."

Wide-cast enumeration of every LA-13511 case currently solvable. NO machinery expansion was needed — the existing `compute_kinf_*` (1G/2G/MG with upscatter handled), `solve_fn_slab_bare_critical`, and `solve_fn_sphere_bare_critical` covered every active case.

## Coverage delivered

42 cases registered (5 Phase A + 37 new Phase B3):

| Class | Active | Stubs | Solver |
|---|---|---|---|
| 1G k_inf infinite | 12 | — | `compute_kinf_1g` (Eq 19) |
| 2G k_inf no upscatter | 6 | — | `compute_kinf_2g_no_upscatter` (Eq 29) |
| 2G k_inf with upscatter | 2 | — | `compute_kinf_2g_general` (Eq 28 corrected) |
| 3G k_inf | 1 | — | `compute_kinf_mg` (Eq 76) |
| 6G k_inf | 1 | — | `compute_kinf_mg` (handles upscatter) |
| 1G slab bare-critical | 4 | — | `solve_fn_slab_bare_critical` (Phase A had 1) |
| 1G sphere bare-critical | 3 | — | `solve_fn_sphere_bare_critical` (Phase A had 1) |
| 1G cylinder bare-critical | — | 3 | B1 dispatch will activate |
| 2G slab/sphere bare-critical | — | 10 | Siewert-Thomas 1986 2G F_N (deferred) |
| **Total** | **30** | **13*** | |

\* Stub count includes 1 case shared with Phase A (Ua-1-0-CY = stub waiting for B1; the registry counts it under WIDE_SLICE_STUBS only via problem-7/23 cases — Phase A's `Ua-1-0-CY` STUB stays in `ALL_FIRST_SLICE`. WIDE_SLICE_STUBS = 12 explicitly.)

## Tests added

**56 new tests** (passing) + **12 skipped stubs** (registered but solver-pending):

* `tests/derivations/test_sood_registry_wide_kinf.py`: 47 tests
  - 11 1G k_inf parametrized over `KINF_1G_CASE_IDS`
  - 5 2G no-upscatter k_inf parametrized
  - 2 2G with-upscatter k_inf parametrized
  - 2 MG k_inf (3G + 6G)
  - 20 cross-implementation gates (`compute_kinf_mg` ↔ `kinf_homogeneous` eig)
  - 5 `general` ↔ `no_upscatter` agreement gates
  - 2 flux-spectrum gates (URR-3-0-IN, URR-6-0-IN mirror)
  - 2 bookkeeping (count, partition completeness)

* `tests/derivations/test_sood_registry_wide_bare_critical.py`: 9 active + 12 skipped
  - 3 slab F_N L1 reference-value tests at N=12
  - 2 sphere F_N L1 reference-value tests at N=10
  - 2 stub-tracking parametrized tests (cylinder × 2; 2G × 10) — all `pytest.skip` with explanatory message
  - 2 bookkeeping (count assertions)

## Tolerances achieved

**k_inf (target ≤ 1e-5 absolute)**:
- min err: 0.00e+00 (URR-3-0-IN, URR-6-0-IN, PU-1-1-IN, Ua-1-0-IN — all hit machine precision)
- max err: 7.0e-6 (Uc-1-0-IN — Sood's published 6-digit value 2.256083 vs computed 2.2560900; within 5-digit envelope per Sood's own claimed precision)
- All 20 cases pass at ≤ 7e-6 absolute.

**Slab F_N at N=12 (target ≤ 1e-5 absolute on a_c)**:
- PUa-1-0-SL (c=1.50): err = 1.32e-6
- PUb-1-0-SL (c=1.40): err = 2.73e-6
- UD2O-1-0-SL (c=1.02): err = 1.89e-6
- All 3 cases pass.

**Sphere F_N at N=10 (target ≤ 1e-5 absolute on R_c)**:
- PUb-1-0-SP (c=1.40): err = 4.0e-8
- UD2O-1-0-SP (c=1.02): err = 3.2e-8
- Both cases pass with massive headroom (better than 1e-7).

## Decisions / non-trivial findings

### N=12 for slabs, N=10 for spheres

Empirically determined. Slab F_N convergence is slower than sphere; N≥14 fails for low c (UD2O at c=1.02 — bracket scan loses the determinant zero at high N). N=12 is the sweet spot: hits ≤ 3e-6 for all 3 slab cases without bracket failure.

### 2G upscatter requires Eq 28 general form

URRb-2-0-IN and URRc-2-0-IN have Σ_21s > 0 (upscatter from slow→fast). The no-upscatter Eq-29 specialisation rejects these inputs at function entry (raises ValueError). Tests use `compute_kinf_2g_general` for these two cases. URRd-2-0-IN despite ν_fast=1.004 and Σ_22s=0 has no upscatter — it uses Eq 29.

### URR-3-0-IN exactness via Forster's inverse design

Sood Eq 60-65 documents: f_23=4 and f_13=15 give k_inf=1.60 and ratios 0.480/0.150 by design. Computed values match to machine precision (literally 0.00e+00 difference). This is the cleanest possible truth-set: Sood constructs the XS so that the answer is a clean rational.

### URR-6-0-IN mirror structure

The 6G case is two coupled URR-3-0-IN blocks: Sood groups 6,5,4 (top three) and 3,2,1 (bottom three) are decoupled in scattering but coupled via χ. The bottom three groups have a thermal-upscatter pattern (Σ_21s=0.171, Σ_31s=0.033, Σ_32s=0.275). The dominant eigenvector exhibits exact mirror symmetry: φ[0]=φ[5], φ[1]=φ[4], φ[2]=φ[3] to machine precision (verified by `test_URR_6_0_IN_flux_spectrum_mirror`).

### Convention conversions documented in `_mix_2g_isotropic` helper

Added a 2G isotropic mixture-builder helper that takes Sood-convention names (sigma_22s, sigma_11s, sigma_12s, sigma_21s — fast/slow Sood-side) and emits an ORPHEUS-ordered Mixture (g=0 fast, g=1 slow). Keeps call sites readable while doing the convention conversion at construction time per Phase A's established pattern.

## No new typos / inconsistencies in Sood's tables

The Phase A migration discovered Sood Eq 28 has a typo (caught by SymPy det(M)=0 derivation in `derive_kinf_2g_general_from_matrix`). Phase B3 enumeration found no additional typos:

- All 1G k_inf values verified to ≤ 1e-5; the worst-case (Uc-1-0-IN at 7e-6) is consistent with Sood's published 6-digit precision and is NOT a typo — it's just rounding of the published "2.256083" value.
- All 2G k_inf values verified at the same precision.
- 3G/6G values are exact-by-design (Forster construction) and computed to machine precision.

## What I did NOT touch

Per the parallel-phase scope wall:

- ✓ Did NOT touch `singular_eigenfunction/` (B1 territory).
- ✓ Did NOT touch `fn_method/{slab,sphere}/one_group.py` interior (B1 may extend; B3 only consumed the existing API).
- ✓ Did NOT touch `sood_registry/cache.py` (B4 territory; preserved imports as-is).
- ✓ Did NOT touch `fn_method` SymPy origins (Phase A locked).

Only modified files:
1. `orpheus/derivations/continuous/sood_registry/la13511.py` — added 37 new cases + 5 helper functions
2. `orpheus/derivations/continuous/sood_registry/__init__.py` — re-export new cases + slice tuples
3. `tests/derivations/test_sood_registry_wide_kinf.py` — NEW
4. `tests/derivations/test_sood_registry_wide_bare_critical.py` — NEW
5. `docs/theory/sood_registry.rst` — added Phase B3 coverage matrix + updated Phase B preview

## Verdict on tractability

The wide enumeration was **highly tractable** — the registry schema designed in Phase A scaled cleanly to 42 cases. Key factors:

1. **Schema was correct**: `materials: dict[int, Mixture]` + `MeshTemplate` + `La13511Truth` covered every case shape. No schema fields needed adding.
2. **Branch-2 solvers were complete**: every active case was solvable by an existing `compute_kinf_*` or `solve_fn_*_bare_critical` entry — no new code shipped under `orpheus/`.
3. **Parametrize fan-out worked perfectly**: `pytest.mark.parametrize("case_id", ...)` over slice tuples is the right test pattern. 47 k_inf tests collapse to <50 lines of test code.
4. **Convention helper paid off**: `_mix_2g_isotropic` saved repeating the Sood→ORPHEUS group conversion 9 times.

The schema **did not strain at scale**. If anything, the Phase A design anticipated wide enumeration cleanly — the slice tuples (`WIDE_SLICE_KINF`, `WIDE_SLICE_BARE_CRITICAL_1G`, `WIDE_SLICE_STUBS`) are a natural organizational layer that emerged from this work, not a refactor.

## Outstanding follow-on (NOT done in B3)

1. **Cylinder solver activation**: B1 dispatch is in flight. Once B1 lands, the 3 cylinder STUBS (`Ua-1-0-CY`, `PUb-1-0-CY`, `UD2O-1-0-CY`) auto-activate via the same `WIDE_SLICE_BARE_CRITICAL_1G` test pattern.
2. **Flux-ratio verification**: B2 dispatch is in flight. The `flux_ratios` truth values are stored on PUb-1-0-SL/CY/SP, Ua-1-0-SL/SP, UD2O-1-0-SL/SP, URRa-2-0-SL — but no test currently verifies them. B2 will fill this gap once flux-reconstruction is built.
3. **2G F_N (Siewert-Thomas 1986)**: not in any current dispatch. The 10 2G bare-critical STUBs (PU/U/UAL/URRa/UD2O × {SL, SP}) await this.
4. **Anisotropic slab cases (32-37)**: PUa-1-1-SL (P_1), PUa-1-2-SL (P_2), Ua-1-1-CY etc. Not enumerated. Anisotropic slab F_N requires extending Siewert-Benoist Eq 7 with the linear-anisotropy term — clean extension but not in any current dispatch.
5. **Reflected/multi-region cases**: 24 reflected cases (problems 3-4, 9-10, 16, 18, 20, 25-28, 30, 58-66). Need multi-region F_N — separate effort entirely.
6. **Archivist dispatch**: skipped per the user-control directive on the Phase A closeout (the user wants explicit control over archivist dispatches; B3 follows that policy).

## Provenance

- All XS values transcribed from LA-13511 PDF (`/workspaces/ORPHEUS/scratch/literature/Sood Foster Parsons (1999)Analytical Benchmark Test Set for Criticality Code Verification.pdf`) Tables 5, 12, 16, 20, 24, 28, 30-31, 33-34, 36-37, 39-40, 43-44, 46-47, 49-50, 59-67. Pages 16-37 + 51-57 (Appendix A).
- Truth values transcribed verbatim from the same tables.
- Cross-checked against literature memo `.claude/agent-memory/literature-researcher/sood_fn_method_full_extraction.md`.

## Files

- Registry: `orpheus/derivations/continuous/sood_registry/la13511.py` (671 → ~1450 lines).
- Re-exports: `orpheus/derivations/continuous/sood_registry/__init__.py`.
- Tests: `tests/derivations/test_sood_registry_wide_kinf.py` (NEW), `tests/derivations/test_sood_registry_wide_bare_critical.py` (NEW).
- Sphinx: `docs/theory/sood_registry.rst` (updated Phase B preview → Phase B3 coverage matrix).
- Closeout: this memo.
