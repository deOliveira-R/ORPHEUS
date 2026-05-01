# `orpheus.derivations`

Reference solutions and symbolic discretisations that anchor the
ORPHEUS V&V ladder. The package is organised into **two verification
paths** that branch from a **shared mathematical root**:

```
derivations/
├── common/        ← math root + cross-cutting utilities
├── discrete/      ← Path 1: symbolic discretisations (production solvers)
└── continuous/    ← Path 2: continuous reference solutions
```

Both paths start from the same canonical first-order transport
equation (laid out symbolically in
`common/transport_equation.py`). Path 1 applies the explicit
discretisation each production solver commits to; Path 2 uses
analytical, semi-analytical, or transcendental techniques that bypass
spatial discretisation. Equality of the two paths at the eigenvalue
or flux-shape level is what verifies the solver.

## Top-level layout

### `common/`

Shared utilities that no single path owns. The four most-imported:

- `kernels.py` — exponential-integral $E_n$, Bickley–Naylor $\mathrm{Ki}_n$,
  chord-half-length primitives.
- `quadrature.py` — unified `Quadrature1D` value object plus the
  primitive constructors (Gauss–Legendre, composite, visibility-cone,
  Gauss–Laguerre) and `AdaptiveQuadrature1D` for `mpmath.quad` rules.
- `quadrature_recipes.py` — geometry-aware recipes
  (`chord_quadrature`, `observer_angular_quadrature`).
- `xs_library.py` — synthetic cross-section library shared between
  derivations and tests.

Plus the type root (`verification_case.py`,
`continuous_reference.py`), the eigenvalue helpers (`eigenvalue.py`),
the shifted-Legendre primitives (`shifted_legendre.py`), and the
symbolic-transport-equation contract module
(`transport_equation.py`).

### `discrete/` — Path 1: symbolic production-solver discretisations

Each sub-package mirrors a production solver in `orpheus/` and
holds the **explicit symbolic form** of the discrete equations the
solver must satisfy.

- `discrete/sn/balance.py`, `discrete/sn/contamination.py` —
  S$_N$ balance + Morel-Montry contamination.
- `discrete/moc/equations.py` — MOC characteristic-step equations.
- `discrete/cp/` — empty placeholder (Issue #141: lift the
  symbolic CP discretisation out of the continuous flat-source CP
  derivations once the path split is consummated).
- `discrete/diffusion/` — empty placeholder; the analytical
  diffusion gap is tracked alongside the CP gap.

### `continuous/` — Path 2: continuous reference solutions

References that are mesh-independent — the production solvers
converge to them under refinement.

- `continuous/analytical/homogeneous.py` — closed-form
  infinite-medium $k_\infty$ for synthetic XS sets.
- `continuous/flat_source_cp/{slab,cylinder,sphere,geometry}.py` —
  flat-source CP eigenvalues (the legacy reference for the
  flat-source production CP solver in `orpheus.cp`).
- `continuous/peierls/` — *the active heterogeneous reference for
  CP verification*. Peierls integral form solved via high-precision
  Nyström collocation at 30+ digits. Geometry-specific solvers
  (`slab.py`, `cylinder.py`, `sphere.py`) plus shared infrastructure
  (`geometry.py`, `cases.py`, `reference.py`). Symbolic *origins*
  (specular-BC R-matrix, cylindrical 3-D G-BC, Knyazev shifted-
  Legendre identities) live under `continuous/peierls/origins/`.
- `continuous/greens_function/` — empty placeholder; planned in
  Plan-2 for fixed-source analytical Green's-function references.
- `continuous/mms/{sn,moc}.py` — Method of Manufactured Solutions
  for differential transport (S$_N$ + MOC), the L1 reference target
  the verification campaign migrated SN/MOC heterogeneous cases to.
- `continuous/cases/{sn,moc,mc,diffusion}.py` — per-method
  verification-case registries that compose the derivations above
  into a name-keyed manifest the top-level
  `reference_values` registry pulls from.

### Top-level files

- `__init__.py` — re-exports the top-level retrieval API
  (`get`, `continuous_get`, …).
- `reference_values.py` — central case registry. Lazy-loaded; both
  the legacy `VerificationCase` and the Phase-0
  `ContinuousReferenceSolution` registries live here. Walks every
  `orpheus.derivations.*` sub-package looking for
  `continuous_cases()` entry points (auto-discovery — no manual
  registration).
- `generate_rst.py` — emits the `docs/_generated/` reference tables
  consumed by Sphinx.

## Notes

### Richardson extrapolation

The Richardson-extrapolation reference path was retired in code by
2026-04. All heterogeneous reference cases now use analytical or
semi-analytical techniques (Peierls Nyström for CP, Method of
Manufactured Solutions for SN/MOC, transfer-matrix for diffusion).
The `common/_richardson_cache.py` utility is retained as a generic
reference cache; it no longer caches Richardson extrapolations
specifically. A future cleanup may rename it to
`common/_reference_cache.py`.

### Discrete-path placeholders

`discrete/cp/` and `discrete/diffusion/` are intentionally empty
sub-packages today. The CP gap is tracked as Issue #141 — the goal
is to lift the symbolic discrete-CP form out of the continuous
flat-source CP derivations, mirroring the SN balance / MOC equations
modules. The diffusion gap is analogous and will be addressed
alongside the CP work.

### Greens-function references

`continuous/greens_function/` is reserved for Plan-2 work. The intent
is to provide closed-form analytical Green's-function references for
fixed-source transport problems (point/line sources in an infinite
medium) — useful as L1 reference targets for the SN and MOC sweep
verification ladders.

## See also

- `docs/theory/verification.rst` — full architectural treatment of
  the two-path V&V architecture and why the `discrete/` ↔
  `continuous/` split exists.
- `docs/verification/reference_solutions.rst` — the contract that
  `ContinuousReferenceSolution` commits to (operator form,
  problem spec, provenance).
- `docs/api/derivations.rst` — Sphinx-generated reference for the
  packages and modules listed above.
