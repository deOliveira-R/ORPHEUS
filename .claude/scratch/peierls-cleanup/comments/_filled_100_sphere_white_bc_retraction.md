# Issue #100 close-out — sphere white-BC rank-1 retraction (historical record)

**Date:** 2026-04-30
**Doc cleanup commit:** [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits)
**Source relocated from:** `docs/theory/peierls_nystrom.rst` lines 3824–3915 (and a duplicate copy at `docs/theory/collision_probability.rst:3922–4051` handled separately).

---

## Summary

- The original Issue #100 conclusion ("rank-1 fails structurally on the sphere because P_esc/G_bc varies ~40 % across radius, giving k_eff ≈ 6.7") was an **artifact** of a missing R² surface divisor (the cylinder rank-1 code path was reused for the sphere without updating `A_d = 2πR → 4πR²`).
- After the divisor fix dispatched by `CurvilinearGeometry.rank1_surface_divisor`, the sphere rank-1 white-BC closure gives **physically sensible behaviour matching the cylinder** (e.g. 0.7 % error at R = 5 MFP, 0.3 % at R = 10 MFP).
- **Production decision (preserved): the rank-1 sphere white-BC closure is structurally correct.** The original "ratio varies 40 %" argument conflated the rank-1 outer product `u_i v_j` with the per-current `J⁻` and is wrong as a structural-failure argument. Sphere Peierls is shipped (Phase-4.3 commits `435c0b3`, `9d03948`, `cad2f0b`); the Peierls-vs-CP comparison at white-BC parity runs in `tests/cp/test_peierls_sphere_flux.py`.

Surviving doc anchors:
- `:ref:` `issue-100-retraction` in `collision_probability.rst` (CP-page version of the retraction debate; the unified-page version is now this comment + a stub).
- `:ref:` `peierls-rank-n-bc-closure`, `peierls-rank-n-jacobian-derivation`, `peierls-rank-n-P-esc-moment` in `peierls_nystrom.rst` (the rank-N skeleton that survives the §8 trim).
- `:label:` `peierls-M-rank-1`, `peierls-M-rank-2` (Lambert/Marshak change-of-basis matrices retained in §F.4 / §F.5).

<details><summary>Full investigation record (verbatim from peierls_nystrom.rst:3824–3915)</summary>

### 2026-04-18 update — retraction

The Phase-4.3 unified sphere Peierls implementation delivers **physically sensible rank-1 white-BC behaviour matching the cylinder**. The "$k_{\rm eff} \approx 6.7$" datum and the "rank-1 fails structurally on the sphere" conclusion below are artefacts of the earlier attempt's missing $R^{2}$ surface divisor (the cylinder code was repurposed for the sphere without updating $A_d = 2\pi R \to 4\pi R^{2}$). The corrected divisor is now dispatched by `CurvilinearGeometry.rank1_surface_divisor`.

Current sphere rank-1 $k_{\rm eff}$ scan (bare sphere, $\Sigma_t = 1$, $\Sigma_s = 0.5$, $\nu\Sigma_f = 0.75$, $k_\infty = 1.5$):

| R/MFP | $k_{\rm eff}$ | err vs $k_\infty$ |
|-------|---------------|-------------------|
| 1.0   | 1.0963        | 26.9 %            |
| 2.0   | 1.3914        | 7.2 %             |
| 5.0   | 1.4897        | 0.7 %             |
| 10.0  | 1.4957        | 0.3 %             |
| 20.0  | 1.4945        | 0.4 %             |

This **parallels the cylinder** (21 % at $R=1$ MFP, falling to 1 % at $R=10$ MFP): both geometries show the same inverse-cell-size growth of the rank-1 Mark closure error, which is a flat-source artefact reapplied at the pointwise level (Issue #103 / N1). The full retraction discussion and the R-vs-R² gotcha are archived in `collision_probability.rst` § `issue-100-retraction`. The text below is preserved as historical context — keeping the record of what was tried and why it failed prevents the same mistake from being made twice.

Sphere Peierls is **shipped** as of Phase-4.3 (commits `435c0b3`, `9d03948`, `cad2f0b`); the Peierls-vs-CP comparison at white-BC parity runs in `tests/cp/test_peierls_sphere_flux.py`.

### Historical text (pre-correction)

The identical failure mode was observed in the Phase-4.3 spherical Peierls attempt (GitHub Issue #100). The sphere's uncollided escape probability $P_{\rm esc}(r)$ varies from ~0.37 at the centre to ~0.68 at the surface, while the re-entry distribution $G_{\rm bc}(r)$ varies from 0 at the centre (Davison's $u(0) = 0$ constraint) to ~2.7 at the surface, and the ratio is not constant — it varies by ~40 % across the sphere radius. A rank-1 correction necessarily imposes a *constant* ratio, so it over-shoots near the surface and under-shoots near the centre, giving $k_{\rm eff} \approx 6.7$ for a 1-G 1-region case (expected $k_\infty = 1.5$).

Both observations — the cylinder's size-dependent error and the sphere's structural failure — are **the same phenomenon**: rank-1 is a flat-source result re-applied at the pointwise level. The two paths forward (and open as of the session of this writing) are:

(a) *Augmented Nyström system*: add the surface-current unknowns as additional degrees of freedom, promoting the $(N\times N)$ system to $(N+n_{\rm surf})\times (N+n_{\rm surf})$. The rank of the white-BC block grows from 1 to $n_{\rm surf}$, which represents the angular resolution of the re-entering distribution.

(b) *Higher-rank angular decomposition*: resolve $J^{-}$ in a Mark-$n$-like $P_n$ expansion of the re-entering hemisphere. Rank $n+1$ correction.

### Post-correction assessment (2026-04-18)

The "ratio varies by 40 %" argument above conflates two independent things: the rank-1 closure is an outer product $u_i\,v_j$ where $u$ and $v$ can individually vary with radius. What the rank-1 closure approximates is the re-entering **angular distribution** $J^{-}(\Omega)$ (treated as uniform isotropic by Mark), **not** the $(i, j)$ coupling structure. A radius-dependent ratio $P_{\rm esc} / G_{\rm bc}$ therefore does **not** imply structural failure; it is absorbed into the outer-product factorisation. What the rank-1 closure actually suffers from is the Mark-closure error in the angular shape of $J^{-}$, and that error scales with cell optical thickness (thick cells homogenise the angular distribution via multiple scattering, thin cells do not). Path (a) and path (b) above are the correct architectural fixes — Issue #103 (N1) tracks higher-rank angular decomposition — but they apply **equally** to cylinder and sphere. Neither is a sphere-specific blocker.

</details>

---

## Why this comment exists

Per the post-#138 documentation cleanup directive, failed-experiment narrative belongs in GitHub issues, not in evergreen Sphinx theory pages. The doc retains a ~25-LoC stub pointing here for the historical debate; the production rank-1 sphere white-BC behaviour and its cylinder-parallel error scan stay in the surviving §8 white-BC closure section (the `peierls-rank-n-jacobian-derivation` label is preserved).

If a future agent picks up the path-(a) augmented-Nyström or path-(b) higher-rank decomposition for white-BC, this comment is the load-bearing record of why neither is a *sphere-specific* blocker — the cylinder shares the same $1/R$ MFP scaling of the Mark closure error.
