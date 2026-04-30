# Archived modules — `orpheus.derivations`

Code removed from the active `orpheus.derivations` package but preserved
verbatim because it may be useful for future work. Each archive entry
records *why* it was removed and *when* it should come back.

## Directory layout

```
derivations/archive/
├── README.md                              ← this file (the index)
├── peierls_moments.py                     ← see "Moment-form Nyström"
├── peierls_slab_moments_assembly.py       ← slab K assembly via moments
├── peierls_cylinder_polar_assembly.py     ← cylinder-polar (φ-quad)
└── peierls_class_b_sphere_bickley_naylor.py ← see "Sphere class-B Bickley-Naylor"
```

The `tests/` subdirectory was removed on 2026-04-30 — both
`test_peierls_moments.py` and `test_peierls_slab_moments.py` had
bit-rotted (orphaned imports of `orpheus.derivations.peierls_moments`
and `slab_polar_K_vol_element` after the post-#131 consolidation in
commit 529cdbe). When the moment-form path is revived for issue #117,
the tests should be regenerated against the lifted module rather than
restored from git history.

## Archive entries

### Moment-form Nyström for production CP — see [GitHub Issue #117](https://github.com/deOliveira-R/ORPHEUS/issues/117)

**Files:**
- `peierls_moments.py` — closed-form polynomial moments for $E_n$,
  $\mathrm{Ki}_n$, $e^{-u}$.
- `peierls_slab_moments_assembly.py` — slab K assembly via the
  moment-form Nyström architecture (closed-form, exact, ~200 ms for
  N=24 at p=6, dps=30).
- `tests/test_peierls_moments.py` — 32 L0 gates of the closed-form
  moment recursions vs `mpmath.quad` (1e-13 / 1e-12 / 1e-15 tolerances
  per kernel family).
- `tests/test_peierls_slab_moments.py` — 9 L1 gates of the slab
  moment-form K matrix vs the legacy E_1 Nyström and the adaptive
  polar reference (1e-12 / 1e-10 elementwise).

**Why archived:** Issue #117 captures the full architecture and
literature derivations. The verification side (`peierls_geometry`) does
not need a fast K assembly — it can afford adaptive `mpmath.quad` per
element (`K_vol_element_adaptive`). The moment form is the *production*
path for a future higher-order discrete CP solver; this branch's CP
production module (`orpheus.cp`) uses flat-source CP and does not need
it.

**When to bring back:** When implementing higher-order spatial source
expansion in production CP, or when adding LS-MOC / quadratic-MOC to a
future MoC production solver. See Issue #117 for the trigger
conditions and the conditioning / Vandermonde caveats.

### Cylinder-polar assembly (explicit out-of-plane φ quadrature)

**Files:**
- `peierls_cylinder_polar_assembly.py` — the
  `_build_volume_kernel_cylinder_phi` body that computes
  $\mathrm{Ki}_1(\tau)$ via $\sum_{\varphi} \exp(-\tau/\cos\varphi)\,w_\varphi$
  with a 16-point GL on $[0, \pi/2]$.

**Why archived:** Mathematically equivalent to `cylinder-1d`
($\mathrm{Ki}_1$ evaluated via `ki_n_mp` / `ki_n_float` directly).
Verified element-wise to machine precision against `cylinder-1d` at
n_phi=32. The cylinder-polar route was a detour in the
"retire-Bickley" sub-thread of issue #116 — useful for the
exposition of the φ-decomposition but not a separate physics
construct.

**When to bring back:** If a future analysis needs to expose the
out-of-plane angular distribution explicitly (rather than as an
integrated Bickley value), e.g., for higher-order angular flux
moments at the cell surface, this is the assembly to start from.
Otherwise, `cylinder-1d` with closed-form `ki_n_mp` is the natural
form.

### Sphere class-B Bickley-Naylor SymPy derivation — see [GitHub Issue #101](https://github.com/deOliveira-R/ORPHEUS/issues/101)

**Files:**
- `peierls_class_b_sphere_bickley_naylor.py` — first-principles
  SymPy derivation of solid-sphere $P_\text{esc}(r)$ via the angular
  $\theta$ form, $u = \cos\theta$ substitution, and $t = \rho$
  recognition that lands at the closed-form $E_n$ expression.

**Why archived:** Stalled investigation. The derivation reaches a
clean closed form, but the production primitive it was meant to feed
(`compute_P_esc_bickley` / `compute_G_bc_bickley`) was never written
in `orpheus.derivations.peierls_geometry`. Without a Sphinx breadcrumb
or a paired test, the file was a Cardinal Rule 1 violation
(orphan derivation at the package root). Archived 2026-04-30 with the
expectation that Issue #101 / #132 will revive it as the math origin
when the sphere class-B closure is closed-out.

**When to bring back:** When implementing
`compute_P_esc_bickley(r)` and `compute_G_bc_bickley(r)` in the
production package. The derivation is the **canonical
specification** — the production code should import from a lifted
`orpheus.derivations.peierls_sphere` and a paired
`tests/derivations/test_peierls_sphere_bickley_symbolic.py` should
verify production-vs-symbolic at machine precision (template:
`tests/derivations/test_peierls_specular_symbolic.py`, established
2026-04-30 as the math-origin → bifurcation pattern POC for
issue #95).
