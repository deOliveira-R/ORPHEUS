---
name: Three meanings of "transport Green's function" — Sanchez/Chandrasekhar gap analysis (2026-05-03)
description: Disambiguates Meaning 1 (Chandrasekhar angular G(τ,τ';μ,μ')), Meaning 2 (Sanchez source-reduction kernel K(r,r')), Meaning 3 (trajectory MoC multi-bounce). Identifies which local PDFs cover which meaning. Local literature is sufficient for slab + sphere with isotropic OR linearly anisotropic scattering. Cylinder Meaning 2 is the only path requiring acquisition.
type: project
---

# Three meanings of "transport Green's function" + literature mapping

When the user asks about Sanchez or Chandrasekhar Green's functions,
disambiguate first:

## Meaning 1 — Chandrasekhar angular Green's function `G(τ,τ';μ,μ')`

Response to a δ-source in space AND direction. Built via singular
eigenfunctions + H-function (or X-function for finite media).
**No ray tracing.**

**Local sources sufficient**:
- Mitsis 1963 ANL-6787 (full PhD thesis, Ch. II-V)
- Atalay 1997 ProgNuclEng 31, 229-252 (anisotropic linear, slab + sphere)
- Kaper-Lindeman-Leaf 1974 NSE 54, 94 (Fredholm `A(ν)`/`B(ν)` recipe)
- Westfall-Metcalf 1973 NSE 52, 1-11 (cylinder)
- Siewert-Benoist 1979 + Grandjean-Siewert 1979 (F_N theory + numerics)

**Already implemented in ORPHEUS** (slab + sphere isotropic):
`fn_method/{slab,sphere}/flux_reconstruction.py`. Gap is API exposure —
expose `_X_negative_real_axis` as public, add fixed-source version of
`solve_kll_*_continuum_coefficient`, and surface `G_chandrasekhar(...)`
as a top-level callable.

**Already implemented for slab + sphere anisotropic**:
`case_method/{slab,sphere}/one_group.py` (Atalay 1997 Eqs. 46/54
criticality). Gap: it's a determinant solver, not a Green's function —
needs a fixed-source extension.

**Cylinder gap is real** (not just API): `singular_eigenfunction/cylinder/`
has WM-72 critical-radius solver but not the interior `ψ(r, μ)`
reconstruction analog of slab/sphere `*_angular_flux_from_scalar`.

## Meaning 2 — Sanchez source-reduction kernel `K(r, r')`

Closed-form spectral kernel: `φ(r) = ∫K(r,r')q(r')dr'` where K is
explicit `E_n` / `Ki_n` / X-function combinations after spectral
integration over μ. **No ray tracing, no eigenfunction expansion at
runtime**.

**Local sources sufficient for sphere**:
- **Sanchez 1986 TTSP 15, 333-343** (sphere with linear anisotropy +
  specular + diffuse BC) — Eqs. (A1)-(A7) appendix is the closed-form
  kernel.
- **Pomraning-Siewert 1982 JQSRT 28, 503-506** (isotropic precursor) —
  Eq. (21) is the explicit Fredholm form `r·I(r) = r·G(r) + ∫₀^R x[(ω/2)I(x) + Q(x)][E_1(|r-x|) - E_1(r+x) + α·F_1(r,x) + β·F_2(r,x)]dx`.
  Eqs. (22)-(25) give G, F_1, F_2 in closed-form 1D quadratures.

**ORPHEUS status**: NOT implemented as a separate spectral evaluator.
Variant α in `peierls_greens_function/` evaluates the SAME kernel
(Sanchez Eq. A5) but via ray-traced multi-bounce, NOT via the closed-form
1D μ-quadrature of Eq. (A6).

**Cylinder Meaning 2 is the only acquisition gap**: Sanchez 1986
covers sphere only (Table 1 lists `α=1` for cylinder but the paper does
not derive it). For cylinder Meaning 2 the user would need:
- Carlvik 1965 Geneva Conf. Vol.2 p.225 (NOT LOCAL — only Carlvik 1967
  finite cylinders/cuboids is local)
- OR Knyazev 1993 *Atomic Energy* 74(5):368-374, DOI 10.1007/BF00844623 (sole author A. P. Knyazev — provides closed-form 1-D quadrature for homogeneous linearly-anisotropic cylinder kernel in Bickley-Naylor functions; cylinder analog of Sanchez 1986 sphere Eq A6). LOCAL as of 2026-05-03. [Earlier drafts cited a non-existent "Knyazev-Selivanov 2014 JCP" — that was a hallucinated phantom citation; CrossRef shows no such paper. Corrected 2026-05-03.]

## Meaning 3 — Trajectory MoC multi-bounce closure

Same numerical answer as Meaning 2, computed by ray tracing + geometric
series resolvent `T = (I - S)⁻¹`. This is what ORPHEUS already
implements in `peierls_greens_function/` (Variant α).

## Verification matrix

(γ) Chandrasekhar is structurally independent of (α) Variant α and (β)
Sanchez kernel — different domain (`ν` spectrum vs `r'` characteristic),
different identity (singular-eigenfunction completeness vs Peierls).

(α) and (β) share the `E_1` building block — their cross-check is
weaker (consistency between two evaluations of the same identity).

The triple match (α)≈(β)≈(γ) is L1-grade verification per
`vv-principles`. Stronger than Sood 2003 numerical benchmark.

## When to cite this

If the user asks:
- "What do I need to implement Chandrasekhar?" → Mitsis 1963 + Atalay
  1997 + KLL 1974, all local. Refactor existing `fn_method/.../flux_reconstruction.py`.
- "What do I need to implement Sanchez 1986?" → Sanchez 1986 + PS-1982,
  both local. Build new `peierls_sanchez_kernel/sphere.py`.
- "Cylinder version?" → STOP. Sanchez 1986 doesn't do cylinder. For
  cylinder Meaning 2 — both required references are now LOCAL: Carlvik 1965 + Knyazev 1993 (DOI 10.1007/BF00844623); the older "Knyazev 2014" reference was a phantom citation that does not exist.

## Notation gotcha

Sanchez 1986 Table 1 uses `α ∈ {0, 1, 2}` for slab/cyl/sphere geometry
index. Sanchez 1986 main body uses `α ∈ [0, 1]` for specular reflectivity.
Same letter, two meanings. ORPHEUS convention is `alpha` = specular
reflectivity. Always rename the geometry index when implementing.

## Gap analysis location

Full deliverable: `/workspaces/ORPHEUS/.claude/scratch/sanchez_chandrasekhar_gap.md`
(written 2026-05-03). 8-section gap analysis with verification matrix
and implementation phases.
