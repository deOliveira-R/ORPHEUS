---
name: Chandrasekhar Meaning (1) is already implemented inside fn_method
description: Surprising finding — ORPHEUS already has the Chandrasekhar singular-eigenfunction angular Green's function machinery for slab + sphere; it lives buried inside `fn_method/{slab,sphere}/flux_reconstruction.py` (KLL 1974 path) but is not exposed as a top-level Green's function callable.
type: project
---

# Chandrasekhar (Meaning 1) is latent in `fn_method/.../flux_reconstruction.py`

When asked "what's missing to implement Chandrasekhar's angular Green's
function in ORPHEUS?", the immediate naive answer is "everything — we
have F_N criticality solvers but not the eigenfunction expansion".

That's WRONG. ORPHEUS already has, for slab + sphere with isotropic
scattering:

- `_X_negative_real_axis(mu, c, u0)` — the Wiener-Hopf X-function for
  finite media. Lives in `fn_method/slab/flux_reconstruction.py`.
  Private (underscore-prefixed) but fully working.
- `solve_kll_slab_continuum_coefficient(...)` — KLL 1974 Eqs. 6-7
  Fredholm iteration for `A(ν)`, the half-range expansion coefficient.
  PUBLIC.
- `solve_kll_sphere_continuum_coefficient(...)` — sphere analog,
  KLL 1974 Eqs. 13-15. PUBLIC.
- `slab_angular_flux_from_scalar(...)` — reconstructs `ψ(z, μ)` at
  any interior `(z, μ)` from the converged `A(ν)`. PUBLIC.
- `sphere_angular_flux_from_scalar(...)` — sphere analog. PUBLIC.

The Atalay 1997 path (`case_method/`) ALSO has the X-function — but
it's the c>1 anisotropic Atalay form (Eq. 40), not the isotropic-c<1
KLL form. There are TWO X-function implementations in the repo:

- `case_method/core/x_function.py::atalay_X_function` — anisotropic c>1
- `fn_method/slab/flux_reconstruction.py::_X_negative_real_axis` — isotropic any c

## Why: gap is API exposure, not math

What's actually missing:

1. A top-level callable `G_chandrasekhar(c, a, tau, tau_prime, mu, mu_prime)`
   that takes a δ-source location and returns the response — i.e. the
   "Chandrasekhar angular Green's function" as users expect it.
2. A fixed-source version of `solve_kll_*_continuum_coefficient` (current
   form assumes bare-critical RHS = 0).
3. Promotion of `_X_negative_real_axis` and the discrete/continuum
   eigenfunctions `φ₀±(μ)`, `φ_ν(μ)` from "implicit inside the
   reconstruction call" to "public callables in `case_eigenfunctions.py`".

Why: When the user (or a future sub-agent) plans Meaning (1)
implementation, they should NOT start from scratch. They should refactor
existing `fn_method/.../flux_reconstruction.py` into a fixed-source
Green's-function exposure path. The math is already verified; the V&V
tests already pass; just expose it.

How to apply: When asked about implementing Chandrasekhar Green's
functions, **always check `fn_method/{slab,sphere}/flux_reconstruction.py`
first**. The KLL 1974 machinery is present. The deliverable is a
refactor, not a new derivation.

The cylinder analog (WM-72 boundary `Φ'(μ)` in
`singular_eigenfunction/cylinder/one_group.py`) does NOT yet have the
interior `ψ(r,μ)` reconstruction — that's a genuine gap, not just an
API gap. Required: complete WM-72 Eqs. 30-31 with non-zero source.

## Sphere = antisymmetric slab parity reduction is also already implemented

`case_method/sphere/one_group.py` uses Atalay 1997 Eqs. 47-54 (sphere =
odd modes of slab). Same parity reduction Sanchez 1986 Eq. 5 uses for
the sphere kernel. So when implementing Meaning (2) sphere, the parity
reduction `ḡ_2(ρ'→ρ) = (ρ'/ρ)·[ḡ_0(ρ'→ρ) - ḡ_0(-ρ'→ρ)]` is already
mathematically present in the codebase, just for a different purpose
(criticality determinant vs. spectral kernel). Don't re-derive.
