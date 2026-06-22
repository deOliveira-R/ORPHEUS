---
name: F_N method sphere — Path 2 (true Siewert-Thomas 1986, not PS-1982 wrapper)
description: Replaced PS-1982-wrapper "sphere F_N" stub with true F_N method per Siewert-Thomas 1986. Geometry-sign abstraction unifies slab + sphere F_N. 46→54 fn_method tests; sphere R_c at 3.6e-8 vs Sood Ua-1-0-SP truth (target 1e-5).
type: project
---

# F_N method sphere — Path 2 closeout (2026-05-02)

`feature/peierls-greens-cylinder` branch, HEAD `667ba60`. 6 atomic
commits.

## Goal

Replace the placeholder PS-1982-wrapper "sphere F_N" with a true F_N
method derivation per Siewert-Thomas 1986. Algebra-of-record
discipline: SymPy from first principles must derive sphere F_N
matrix entries AND critical condition AND match Sood `Ua-1-0-SP`
r_c = 2.4248249802 mfp (c = 1.30) at ≤ 1e-5.

## Files

### Created

- `orpheus/derivations/continuous/fn_method/core/fn_matrix.py` (new
  shared assembler).
- `orpheus/derivations/continuous/fn_method/origins/fn_sphere_derivations.py`
  (Branch-1 SymPy, 5 derive_*() functions).
- `orpheus/derivations/continuous/fn_method/sphere/one_group.py`
  (REWRITTEN — was PS-1982 wrapper, now true F_N).
- `tests/derivations/test_fn_la13511_sphere.py` (REWRITTEN — was 4
  SymPy gates + 1 slow Sood gate; now 11 foundation tests).
- `tests/derivations/test_fn_la13511_sphere_xverif.py` (REWRITTEN —
  was 1 PS-1982/Variant α weak cross-check; now 3 L1 cross-checks).

### Modified

- `orpheus/derivations/continuous/fn_method/core/__init__.py`
  (export new assembler).
- `orpheus/derivations/continuous/fn_method/origins/__init__.py`
  (export new sphere derivations, remove old peierls_sphere ones).
- `orpheus/derivations/continuous/fn_method/sphere/__init__.py`
  (point at new solve_fn_sphere_bare_critical).
- `orpheus/derivations/continuous/fn_method/slab/one_group.py`
  (refactored to call shared assembler with geometry_sign=+1).
- `orpheus/derivations/continuous/fn_method/slab/__init__.py`
  (refresh placeholder docstring).
- `docs/theory/fn_method.rst` (replace pivot narrative + 5 new TODO
  stubs for V_fn-sphere-fn.* labels).

### Deleted

- `orpheus/derivations/continuous/fn_method/origins/peierls_sphere_derivations.py`
  (the OLD pivot SymPy module — replaced by fn_sphere_derivations.py).

## Commits

```
3de22cb refactor(fn_method): extract slab/sphere F_N matrix assembler to core/
e37c0d1 feat(fn_method): Branch-1 SymPy derivations for sphere F_N (Siewert-Thomas)
0a290ae feat(fn_method): Branch-2 sphere F_N solver (Siewert-Thomas 1986)
a1e4e53 test(fn_method): foundation + L1 sphere F_N (true F_N, not PS-1982 wrapper)
ce41969 docs(fn_method): replace PS-1982-wrapper sphere narrative with true F_N
667ba60 docs(fn_method): refresh slab/__init__.py docstring after F_N implementation
```

## Test counts before / after

- **Before**: 46 fn_method tests in 89 s (the slow PS-1982 sphere
  gate took ~80 s of bisection-on-R per evaluation).
- **After**: 54 fn_method tests in 11 s. 8 net new tests; runtime
  cut by 8x because sphere F_N solver runs in <1 s vs PS-1982
  wrapper's ~80 s.

Slab tests: 23/23 still pass at the same 4.82e-6 tolerance vs Sood
`Ua-1-0-SL`. The shared-assembler refactor preserved slab behavior
bit-for-bit.

## Achieved sphere R_c tolerance vs Sood `Ua-1-0-SP` (truth: 2.4248249802 mfp at c = 1.30)

| N  | R_c achieved        | Error vs truth | Notes                                       |
|----|---------------------|----------------|---------------------------------------------|
| 5  | 2.4248247664        | 2.1e-7         | Already 100x better than 1e-5 target        |
| 8  | 2.4248248870        | 9.3e-8         |                                             |
| 10 | 2.4248249443        | 3.6e-8         | **Test default**                            |
| 12 | 2.4248249653        | 1.5e-8         |                                             |
| 15 | 2.4248249747        | 5.5e-9         |                                             |
| 20 | 2.4248249719        | 8.3e-9         | Brent precision floor                       |

The default N=10 reaches 3.6e-8 absolute, three decades better than
the 1e-5 target.

## Cross-check tolerance vs Variant α sphere

`test_fn_sphere_vs_variant_alpha_sphere_at_sood_ua_1_0_sp`: F_N
sphere at N=10 predicts R_c = 2.4248249443 mfp; Variant α sphere
at this radius gives k_eff = 1 to within 4.2e-6 (target 1e-5).

This is a **structurally stronger cross-check** than the previous
PS-1982/Variant α comparison: F_N method works in the Case singular-
eigenfunction representation; Variant α works in the bouncing-
characteristic representation. Above the trusted-library line, no
shared in-house code. By contrast, PS-1982 and Variant α both
reduce the Peierls integral equation by different algebraic paths
(procedurally independent only).

## Honest verdict on Branch-1 SymPy at sphere F_N complexity

### Did sphere SymPy derivation hit any new walls beyond what slab had?

**No.** All 5 derivations close in pure SymPy (State 1A —
closed-form). The structural insight — that sphere = slab with one
sign flip — makes the SymPy work mechanical. No `simplify()` hangs,
no eigenvalue-solving above degree 2, no Piecewise compounding. The
heaviest derivation is `derive_sphere_critical_condition` which
expands a 2×2 symbolic determinant — finishes in <1 s.

The literature-researcher's prediction was correct: sphere F_N is
the EASIEST geometry to extend (no Bickley-Naylor / Bessel-K / X
function gymnastics like cylinder F_N would need). Slab and sphere
share the same function basis (shifted Legendre on [0, 1] with
logarithmic kernels).

### Did the X-function geometry-independence verify symbolically?

**Yes.** `derive_x_function_geometry_independence` symbolically
inspects the free symbols of Λ(z) = 1 - cz·atanh(1/z) and the
X-function integrand arg Λ⁺(τ)/(τ-z), confirming that no geometry
parameter (R, a) appears anywhere. This justifies reusing the slab
X-function machinery for sphere verbatim.

The verification has a small implementation note: SymPy's free
symbol introspection compares NAMES, not values. The check
`free_names & geometry_names == set()` is a syntactic-not-semantic
verification, but for the purpose of "the integrand has no R/a
dependence" it's sufficient because the integrand was constructed
from c, z, τ, μ symbols only.

### Was the sign-flip extension as clean as the literature-researcher predicted?

**Yes — cleaner than expected.** The geometry_sign abstraction
condenses the entire slab → sphere extension to a single integer
parameter passed to one assembler function. The Branch-2 sphere
solver is essentially a parallel of the slab solver with two
substitutions:

1. `geometry_sign = -1` (vs slab's +1).
2. **Different collocation grid** — Siewert-Thomas Eq. 38a
   Chebyshev nodes strictly inside (0, 1), vs slab's Grandjean-
   Siewert grid that includes ξ=0 and ξ=1.

The collocation-grid change was NOT predicted in the literature-
researcher memo (the memo cited Eq. 38a but I initially tried to
reuse the slab grid). I learned the hard way that sphere F_N at
ξ=0 has a structural rank deficiency (the geometry-sign-bearing
attenuation row collapses to a constant independent of R AND of s).
The Chebyshev-interior grid avoids this.

## New typos / errata caught

**No new Siewert-Thomas typos found.** The published equations
(Eqs 4, 39, 46) are self-consistent with my SymPy reformulation.
Unlike the Sood Eq 28 finding from Slice 1 (where SymPy det(M) = 0
gave 2.683767 vs the printed 2.862), Siewert-Thomas's algebra
checks out cleanly.

The closest "near-miss" was the Y_α / X_α factor of c convention:
the published Y_α(ν) = c·F_α(ν) carries a leading c, while ORPHEUS
A_α does not. This is NOT a typo — it's a normalisation convention
that doesn't affect det M = 0 (the constant factors out of the
homogeneous system). I documented it in the
`derive_sphere_fn_matrix_entry` docstring so future readers don't
get confused.

## Architectural seams identified for next phase

### Cylinder F_N (Westfall-Metcalf 1972) — separate folder

The next bare-critical geometry on the LA-13511 ramp is the cylinder
(Sood `Ua-1-0-CY`). Per the literature-researcher memo and
Siewert-Thomas's framing, cylinder F_N is **mathematically distinct**
from slab/sphere F_N:

- Slab + sphere F_N: shifted Legendre basis on [0, 1], logarithmic
  kernels, B/A moment recursion.
- Cylinder F_N (Westfall-Metcalf): shifted Bessel basis, Bickley-
  Naylor Ki_n series, modified Bessel K_0 / K_1 special functions.

Cylinder F_N should NOT live in the same folder as slab + sphere
F_N because:

1. The **shared core** (`core/dispersion.py`, `core/moments.py`,
   `core/fn_matrix.py`) is slab/sphere-specific; the moment integrals
   and X-function are based on the Λ(z) = 1 - cz·atanh(1/z)
   dispersion. Cylinder uses a DIFFERENT dispersion function.
2. The **collocation grid prescription** is different (Westfall-
   Metcalf use a different rule).
3. The **geometry_sign abstraction** is slab/sphere-specific —
   cylinder F_N has no such single-parameter generalisation.

**Recommendation for next phase**: cylinder F_N goes in
`orpheus/derivations/continuous/fn_method/cylinder/one_group.py`
with its own dispersion module, moment library, and assembler.
Reuse only the **interface convention** (`SolverResult` dataclass
shape) and the **F_N collocation pattern** (root-finding on det M
= 0). All the math is fresh.

### F_N flux reconstruction — deferred

The current implementation returns the F_N coefficients
`a_α` (eigenvector of M(R_c) at R = R_c) but does not reconstruct
the angular flux ψ(±R, ∓μ) = Σ a_α μ^α nor the scalar flux φ(r).
Sood Case 4 (Ua-1-0-SP) does NOT publish flux ratios, so the test
gates don't need flux reconstruction. But the sphere flux-shape
verification will need it — that's a follow-up slice.

### Sphere F_N at very high N — Brent precision floor

At N ≥ 25, the |det M(R)| at the root drops below 10^-150 (machine
representation limit for `np.linalg.det`). The Brent refinement
hits its precision floor at ~10^-9 absolute on R, giving N=20 the
best practical accuracy. For higher precision, use mpmath-flavored
LU determinant (the Mitsis 1963 / Westfall-Metcalf approach), which
this implementation does NOT do. The Sood 1e-5 target doesn't need
mpmath, so this is documented as future work.

## Sphinx -W build

Clean. 8 new test gates discovered by the V&V audit harness; 1592
total tests collected (was 1555). 2 new equation labels indexed:

- `fn-method-V-fn-sphere-fn-1` through `fn-method-V-fn-sphere-fn-5`
  (5 sphere F_N labels).
- The 4 old `fn-method-V-fn-sphere-1` through `-4` labels were
  removed cleanly when the obsolete `peierls_sphere_derivations.py`
  was deleted.

## DISPATCH_REQUEST emitted

None. The Sphinx page was updated to "stub" quality (TODO markers
on each label per the algebra-of-record discipline). The archivist
agent will be dispatched in a subsequent phase to expand the rich
narrative.

## Lessons for the skill library

### Sphere F_N collocation grid — the rank-deficiency at ξ=0

Documenting this in the `algebra-of-record` skill as a State-1A
edge case. The slab F_N grid (which includes ξ=0) WORKS for slab
because at ξ=0 the slab matrix entry is `B_α(0) + A_α(0)` (the
exp(-∞)·A vanishes, leaving only the B_α(0) row that's actively
constrained by the slab's dispersion-root structure). For sphere
with sign flip, at ξ=0 the entry is `B_α(0) - 0·A_α(0) = B_α(0)`
— SAME row as slab. So the difference between geometries vanishes
at ξ=0. This collapses the rank, masking the sphere root.

Lesson: **reuse a numerical grid only when the math justifies it,
not because the grid is convenient**. The Siewert-Thomas Eq. 38a
Chebyshev-interior grid is the literature's prescription for a
reason. (Trusted-library-line violation: I shared `collocation_points`
with sphere when the math required a different grid.)

### Sphere F_N root detection — |det| minimum, not Re(det) sign change

For supercritical slab (c>1), the imaginary ν_0 makes the matrix
complex but the BC's symmetric structure makes Re(det) cross zero
cleanly at the root. For sphere (anti-symmetric BC), this Hermitian
property is broken — Re(det) does not cross zero clean. Robust
strategy: minimize |det| (or smin) directly via Brent, anchored
by a prominence-filtered first local minimum from a coarse scan.

Lesson: **for problems with broken Hermitian structure, switch
from sign-change to magnitude-min root-finding.** Adding this to
the `numerical-bug-signatures` recognition catalog as a sphere
F_N-specific pattern (or generic to "F_N with anti-symmetric BC").

## Verification coverage achieved

The new sphere F_N has:

- **5 foundation-level SymPy verifications** (algebra-of-record).
- **3 L1 cross-check tests** (vs Variant α, vs Sood truth).
- **1 collocation-grid sanity test** (the Chebyshev-interior
  prescription).
- **1 geometry-sign discrimination test** (smin gap > 3 decades).
- **1 convergence-with-N test** (algebraic convergence verified).
- **1 dispersion-consistency test** (Λ(u_0) = 0 at machine precision).
- **1 coefficient-normalisation test** (a_0 = 1).
- **1 reference-value gate** (Sood truth at ≤ 1e-5).

Total: 11 new tests + 3 L1 cross-checks = 14 new tests for sphere
F_N. The slab F_N retains 23 unchanged tests.
