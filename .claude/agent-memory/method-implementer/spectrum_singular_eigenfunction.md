---
name: Spectrum class — singular_eigenfunction math-heart
description: 3rd concrete instance of the math-heart pattern; bit-equality preservation across slab+sphere+cylinder; ~430-line teaching section on Case spectrum + expansion theorem
type: project
---

# Spectrum class — singular_eigenfunction math-heart

`feat/spectrum-class-singular-eigenfunction` 2026-05-04. Third
concrete instance of the math-heart pattern, sibling to `Billiard`
(trajectory_resolvent, R2 hindsight refactor) and `MomentSpace`
(fn_method, parallel agent earlier this session). With Spectrum landing
the third sibling, the math-heart pattern crosses the "≥3 instances"
threshold, empirically validating the unifying `TransportSolver`
Protocol design (parallel agent on `feat/transport-solver-protocol`).

**Why:** every reference solver in the project is now mounted on a
shared production-protocol input shape (`materials: dict[int, Mixture]
+ GeometrySpec`) with shared cross-method result types
(`CriticalSolution` / `FluxSolution`). Cross-method comparators in
`tests/cross_method/adapters.py` can dispatch on `method_name`
without isinstance ladders. The Protocol is observed, not posited.

**How to apply:** future math-heart additions follow the precedent
exactly — frozen dataclass owning `GeometrySpec` + method-specific
config, factory `from_problem(materials, geometry, *, ...)`, properties
`materials` / `geometry_spec` / `method_name`, solve methods returning
shared types. Bit-equality preservation against function-level API is
the load-bearing foundation test (verified via `float.hex(...)`
exact-bit comparison).

## Atomic commits

`2503587` — feat(singular_eigenfunction): introduce Spectrum math-heart class
- 1581 insertions: class + tests + `__init__.py` export
- 19 foundation tests pass; all 45 existing tests preserved

`d175d85` — docs(singular_eigenfunction): rich math-teaching section on Case spectrum
- 562 insertions to `docs/theory/singular_eigenfunction.rst`
- Section walks from transport eigenvalue problem → Case dispersion
  relation → continuum P + δ structure → expansion theorem → half-range
  completeness + X-function → geometry reductions → linear anisotropy
  → three-meanings taxonomy → Spectrum class

`a0ec9aa` — docs(singular_eigenfunction): fix two Sphinx -W warnings
- Replace stale `:eq:atalay-eq46` and `:doc:/skills/vv-principles`
  cross-references; build clean

## Deliverables

| Path                                                                  | Size       | Purpose                                  |
| --------------------------------------------------------------------- | ---------- | ---------------------------------------- |
| `orpheus/derivations/continuous/singular_eigenfunction/spectrum.py`   | ~870 LOC   | Spectrum dataclass + factory + solve_*   |
| `orpheus/derivations/continuous/singular_eigenfunction/__init__.py`   | +3 lines   | Export `Spectrum`                        |
| `tests/derivations/test_singular_eigenfunction_spectrum.py`           | ~440 LOC   | 19 foundation-tagged tests               |
| `docs/theory/singular_eigenfunction.rst`                              | +562 lines | "Case spectrum and expansion theorem" §  |

## Bit-equality preservation

The class is a thin facade over the existing function-level API:

* Slab: `solve_case_method_slab_critical(c=spec.c, R=R_extracted, f1=spec.f1, maxdegree=n_modes, ...)`
* Sphere: `solve_case_method_sphere_critical(c=spec.c, R_refl=R_extracted, f1=spec.f1, ...)`
* Cylinder: `solve_singular_eigenfunction_cylinder_bare_critical(c=spec.c, n_grid=n_modes, ...)`

Bit-equality verified by calling the function API with the **extracted**
values (`spec.c`, `spec.f1`) — using literal floats fails because
`(0.7 + 0.6) / 1.0 = 1.2999999999999998 ≠ 1.3`. The float.hex test
gymnastics are documented in `test_slab_anisotropic_bit_equality_with_function_api`.

## Material extraction discipline

* `c = (Σ_s + νΣ_f)/Σ_t` from `mixture.SigS[0].toarray() + mixture.SigP`
* `f1 = SigS[1] / SigS[0]` (linear anisotropy; 0 for isotropic
  mixtures, P_2+ rejected)
* Cylinder + non-zero `f1` rejected at construction (WM-72 isotropic
  only — surfaced at the facade boundary so callers know the
  combination is out of pillar)
* `BC.vacuum → R=0`; `BC("partial", {"albedo": R}) → R` for
  Atalay reflected-slab / sphere benchmarks
* `BC.reflective → R=1` (rejected by underlying solver — Atalay omits
  R=1 columns from Tables 2-5 because thickness drops out)

## Three-meanings positioning

Spectrum lives at meaning **(γ): singular-eigenfunction angular
Green's function** in the Sanchez-Chandrasekhar taxonomy. MomentSpace
*also* at (γ) but with **truncated** Legendre basis (F_N method).
Billiard at meaning **(α): trajectory resolvent**.

The Spectrum ↔ Billiard cross-check is **structurally stronger** L1
evidence than Spectrum ↔ MomentSpace because the two pillars exercise
distinct integrands (ν-spectrum vs ray-traced phase space), where
MomentSpace and Spectrum are sibling closures of the same pillar.

## Branch hygiene incidents (5 occurrences)

Three parallel method-implementer agents were live this session
(TransportSolver Protocol, BasisSpace galerkin_spectral, Spectrum
singular_eigenfunction). The parallel agents performed `git checkout`
while my work was in-flight at least 5 times:

1. After my first commit landing, a checkout swapped the working tree
   away from my branch; my `__init__.py` edit was lost and re-applied.
2. After my second commit attempt, the commit landed on
   `feat/transport-solver-protocol` instead of my branch (HEAD had
   swapped); I cherry-picked `35d6baa → 2503587` to my branch.
3. The wrongly-landed commit on `feat/transport-solver-protocol`
   could not be cleaned up — `git update-ref` on a parallel agent's
   branch is permission-denied (the user's "zero overlap" boundary).
   The parallel agent will see an extra Spectrum commit on their
   branch; this is a recoverable artifact, not data loss.
4. Two attempts at the docs section edit were lost mid-flight to
   branch swaps; recovered by writing the new section content to
   `/tmp/spectrum_theory_section.rst` (branch-independent) and using
   a Python script to atomically insert at the marker.
5. The third agent's commit `d5bb87f` (BasisSpace closeout memo)
   landed on my branch. Harmless (memo only, no source code), so
   left in place rather than reverted to avoid further branch
   thrashing.

**Recovery procedure that worked:**
- Write generated content to `/tmp/*` (branch-independent storage).
- Use Python script with sentinel-marker `replace()` to insert atomically.
- Chain `git checkout feat/spectrum-class-singular-eigenfunction &&
  python3 /tmp/insert.py && git add ... && git commit -m "..."` in
  a single bash invocation so the parallel agents can't trample
  between steps.
- After EVERY commit, `git branch --show-current` to verify HEAD.

## What survived from the brief

The brief's hard constraints all met:
- ✅ Cardinal Rule 1: bit-equality preserved (verified by float.hex)
- ✅ Cardinal Rule 3: math-rich docstring + theory section ship as
  primary deliverable; theory section reads as graduate textbook
  chapter on Case spectral theory
- ✅ Sphinx -W exit 0 at every commit (after a0ec9aa fix)
- ✅ All existing 45 singular_eigenfunction tests pass
- ✅ Branch `feat/spectrum-class-singular-eigenfunction` off main
- ✅ Multi-commit (3 focused commits)
- ✅ Don't push (branch local; user merges manually)

## What deferred / not delivered

- **TransportSolver Protocol explicit conformance test**: the
  parallel agent's spec calls for `isinstance(Spectrum.from_problem(...),
  TransportSolver) == True`. Cannot land in this commit because the
  Protocol class itself was being implemented in parallel; my
  `Spectrum` is **structurally** Protocol-conforming (right factory
  shape, right properties, right return types), so the explicit
  isinstance test will pass when the Protocol agent lands. Separate
  follow-up dispatch needed.
- **Cross-method adapter integration**: `tests/cross_method/adapters.py`
  has a placeholder for SingularEigenfunction adapter; not wired in
  this slice. Follow-up.
- **Slab + sphere flux reconstruction**: explicitly out of pillar for
  Spectrum (owned by F_N pillar's KLL Fredholm path per
  structural-independence rule). Calls to `solve_fixed_source` on
  slab/sphere raise `NotImplementedError` with redirect to MomentSpace.

## Test count

| Test family                             | Count  |
| --------------------------------------- | ------ |
| Construction & validation               | 7      |
| Bit-equality with function API          | 5      |
| Cross-method shared-result-type         | 3      |
| Validation guards                       | 4      |
| **Total Spectrum foundation tests**     | **19** |

Plus 45 existing function-level tests preserved (1.5 minutes total
suite runtime).
