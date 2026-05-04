# Folder Naming Taxonomy — Method-Canonical Naming for ORPHEUS

**Date**: 2026-05-03
**Author**: literature-researcher
**Inputs read** (all LOCAL in `/workspaces/ORPHEUS/scratch/literature/`):
Sanchez 2002 (28 pp, full read), Sanchez 1986 (10 pp, full Appendix
read), Pomraning-Siewert 1982 (full 4-page note), Westfall-Metcalf
1973 (full 11 pp), Atalay 1997 (first 10 pp), Mitsis 1963 ANL-6787
(TOC + abstract + first chapter). No paywall hits — all evidence
sourced from the user's own corpus.

---

## Summary

**Folders to RENAME** (4):

| Current | Proposed | Why |
|---------|----------|-----|
| `peierls_greens_function/` | `trajectory_resolvent/` | Sanchez 2002 names the family "**trajectory-based deterministic transport methods**" in its title; the algebraic core is the "scattering resolvent" T = (I−S)⁻¹. "Greens function" is what the method *produces*, not what the method *is*. |
| `case_method/` | `singular_eigenfunction/` (consolidate into existing folder; see Q2) | Author-named — violates rule. Method-canonical name is "**singular eigenfunction expansion**" per Case 1960 / Atalay 1997 / Mitsis 1963 / Westfall-Metcalf 1973 (all four say so verbatim). |
| `carlvik_galerkin/` | `galerkin_spectral/` | "Carlvik" is an author. The method is "**Legendre-Galerkin spectral expansion**" — Carlvik (1968) introduced *recurrences* for matrix elements; Dahl-Sjostrand (1979) introduced the Galerkin matrix construction. Folder should reflect the structural method, not the recurrence library. |
| `sood_registry/` | (KEEP as-is) | Registry of *cases*, not a method — explicitly allowed under user's rule "Author names → registry of cases (acceptable)". |

**Folders to CONSOLIDATE** (1 pair, see Q2):
`case_method/` (Atalay 1997 slab + sphere with linear anisotropy)
+ `singular_eigenfunction/cylinder/` (Westfall-Metcalf 1973
cylinder, isotropic, Bessel-K kernels)
→ unified under `singular_eigenfunction/` with subfolders by geometry.

**Folders to KEEP** (8):

| Folder | Justification |
|--------|---------------|
| `peierls_nystrom/` | Method-canonical: Peierls integral form + Nyström quadrature. Two technical names, both refer to the *method*, not authors. |
| `fn_method/` | Method-canonical: F_N method (Siewert et al.) — established literature shorthand for "moment projection of singular eigenfunctions". |
| `flat_source_cp/` | Method-canonical: collision probability with flat-source approximation. |
| `mms/` | Method-canonical: Method of Manufactured Solutions. |
| `analytical/` | Acceptable: closed-form analytical references (homogeneous medium k_∞ etc.). Could sharpen to `closed_form/` but `analytical` is the standard verification-pillar terminology. |
| `cases/` | Per-method case registries (composes references). Generic but honest. |
| `sood_registry/` | Author-named *case* registry — allowed by rule. |
| `singular_eigenfunction/` | Method-canonical (after consolidation per Q2). |

**New EMPTY folders to RESERVE** (Q4): see table at end.

---

## Q1 — `peierls_greens_function/` rename

### Is the package exactly Sanchez 2002?

**No** — the construction is Sanchez **1986** (Appendix), not Sanchez
2002. The package's own docstring confirms this:

```python
# /workspaces/ORPHEUS/orpheus/derivations/continuous/peierls_greens_function/greens_function.py:88
References
----------
- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
  DOI: 10.1080/00411458608210456.
```

Sanchez 2002 is about *boundary conditions in trajectory-based
methods* — its Eq. (15) is structurally identical to ORPHEUS's
resolvent (see verbatim quote below) but the paper's *focus* is the
exact-vs-approximated treatment of geometrical motions (translation,
rotation, axial symmetry) for 2-D unstructured CP+MoC tracking.

**The kernel ORPHEUS implements is Sanchez 1986 Eq. (A4)** — verbatim,
from the appendix I just read:

> "*I_{in}(µ) = (1/2πA) e^{−τ_{+}} T(µ_{+}) [α δ(µ−µ_{+})/µ_{+} +
>  β/(1−βχ_*) χ(µ) T(µ)]*  (A4)
> *where T(µ) = 1/[1 − α e^{−τ(µ)}]*"

That is the **resolvent T = 1/(1 − α e^{−τ})** the ORPHEUS code uses,
quoted verbatim from a paper that is not author-named but
geometry+content-named (*"Integral form of the equation of transfer
for a homogeneous sphere with linearly anisotropic scattering"*).

### Does Sanchez 2002 prescribe the bouncing-trajectory + multi-bounce closure?

**Yes — explicitly.** Sanchez 2002 Eq. (14)-(15) gives the closure
verbatim:

> "*ψ(x) = ψ_q(L)/(1 − ψ_bd(L)) ψ_bd(x) + ψ_q(x)*  (15)
> *This formula shows how to calculate the angular flux along a
> periodic compound trajectory of length L. It requires the
> simultaneous computation of two angular fluxes: (a) the angular
> flux produced by a unit incoming angular flux in the absence of
> volumetric sources, ψ_bd(x), and (b) the angular flux produced by
> the volumetric sources with zero incoming angular flux, ψ_q(x).*"
> — Sanchez, Mao & Santandrea (2002), Nucl. Sci. Eng. 140, p. 31

The denominator `1/(1 − ψ_bd(L))` IS the rank-1 resolvent T applied
to a periodic compound trajectory, where `ψ_bd(L) = α · exp(−τ_period)`
in the Variant α homogeneous-medium specialisation. ψ_bd is the
boundary kernel; in ORPHEUS notation ψ_q = `B` and the leading
ψ_bd(x) factor is the first-leg attenuation `exp(−τ_first_leg) · α`.

### What does Sanchez 2002 call this method?

The paper title gives the answer: "**Trajectory-Based Deterministic
Transport Methods**". Section III is titled "EXACT TREATMENT FOR
CLOSED DOMAINS" and reads (Sanchez 2002 p. 31):

> "*A periodic trajectory is characterized by the fact that there is
> a one-parameter subgroup G_p ⊆ G that maps the trajectory onto
> itself. As we will see, this fact implies that periodic trajectories
> exist only for angles φ in a subset P of (0, π).*"

The construction Sanchez 2002 uses is what the paper calls **periodic
compound trajectories** with closure via Poincaré recurrence.
ORPHEUS's Variant α is the 1-D radial specialisation: every chord
in a sphere or infinite cylinder is automatically periodic by axial
symmetry — the "compound trajectory" reduces to a single bouncing
chord.

The earlier session's characterisation —

> "Variant α gives the SAME numerical answer as Sanchez's Eq A1
> spectral kernel, but the CONSTRUCTION is structurally MoC: trace
> rays, integrate exp(−Στ), sum bounces via resolvent T = (I−S)⁻¹"

— **matches Sanchez 2002's vocabulary exactly**: trajectory tracking
+ periodic closure via the geometric series in Eq. (15). It just
isn't *Sanchez 2002 specifically* that introduces the construction —
the construction has been in the field since Mitsis 1963 (Smith 1969
PhD thesis is what Westfall-Metcalf cites for exactly the same
periodic-trajectory closure idea), Pomraning-Siewert 1982 derived
it for the sphere with diffuse+specular reflection (`T(µ) = [1 − α
e^{−2Rµ}]^{−1}`, P-S Eq. (14), VERBATIM the same form), and Sanchez
1986 Appendix put it in the modern multi-region notation that the
ORPHEUS code quotes.

### Name proposal — priority list

1. **`trajectory_resolvent/`** *(strongly preferred)* — captures both the
   structural construction (trajectory tracking) and the algebraic
   close-out (resolvent T = (I−S)⁻¹). Aligns with Sanchez 2002's
   own family name "trajectory-based deterministic transport methods"
   and with the cross-domain-attacker's frame name
   "variant_alpha_2surface_bie_frame" (the boundary-integral-equation
   resolvent on 2-surface geometries). The shared core function is
   already named `compute_resolvent_T` and `apply_variant_alpha_closure`,
   so the folder name harmonises with the API surface.

2. `bouncing_characteristic/` — alternative; emphasises the
   characteristic-method viewpoint over the operator viewpoint.
   Defensible but loses the "resolvent" structural content that
   makes the rank-1 / rank-2 distinction visible in the API.

3. `peierls_resolvent/` — keeps the "Peierls integral form ancestry"
   thread alive (ORPHEUS already documents the Peierls heritage in
   the package docstring). But "Peierls" doesn't appear in Sanchez
   1986 / 2002 / PS-1982 — it is the ORPHEUS author's framing, not
   the literature's. Weaker than option 1.

4. `variant_alpha/` — would match the package's own vocabulary
   (`variant_alpha_core`, `apply_variant_alpha_closure`) but
   "Variant α" is **internal ORPHEUS jargon** — no external paper
   uses it. Strongly **disfavoured** under the method-canonical rule.

**My recommendation**: rename `peierls_greens_function/` →
`trajectory_resolvent/`. Fold the existing `variant_alpha_core.py`
into the new folder name as `core.py` (or keep the file name; the
folder name is what the user reads structurally).

### Renaming impact note (Q1)

Scope of touched code:
- **62 Python files** import from the old path
  (`grep -rln peierls_greens_function`).
- **24 test files** in `tests/derivations/test_peierls_greens_function_*.py`
  (rename test file prefixes for consistency).
- **1 Sphinx page** `docs/theory/peierls_greens.rst` references the
  module path 33 times — rename to `docs/theory/trajectory_resolvent.rst`.
- **1 cross-reference** in `docs/theory/peierls.rst` (5 hits) —
  parent page that introduces the family.
- **1 cross-reference** in `docs/theory/singular_eigenfunction.rst`
  (1 hit, just a verification cross-ref).
- **2 Sphinx labels** to update: `_theory-peierls-greens:` →
  `_theory-trajectory-resolvent:` and the `peierls-greens-family-overview:`
  cross-ref label.
- `index.rst` toctree entry: `theory/peierls_greens` →
  `theory/trajectory_resolvent`.

Mechanical refactor — all changes are mechanical replacements
(`peierls_greens_function` → `trajectory_resolvent` and
`peierls_greens` → `trajectory_resolvent` in RST). No semantic
changes to the algebra, kernels, tests, or numerical thresholds.
The `variant_alpha` API names within the module can keep their
internal vocabulary (they're documentation labels, not folder
names) or be renamed for consistency — that's a separate decision
the user should weigh independently.

---

## Q2 — `singular_eigenfunction/` vs `case_method/`

### Q2.a — Same method-class or distinct?

**Identical method-class.** Both are instances of the **Case
singular eigenfunction expansion** (Case 1960, *Ann. Phys. NY* 9, 1).
Direct evidence from the four primary sources (all locally
available, all read):

**Westfall-Metcalf 1973, p. 1, Introduction**:
> "*Since the introduction of the **singular eigenfunction expansion
> technique** by Case [Ref. 1] in 1960, a wide variety of transport
> problems have been treated by this method.*"

**Atalay 1997, abstract**:
> "*Case's singular eigenfunction method is used to formulate the
> criticality conditions. In addition to available bi-orthogonality
> relations in the literature, some parallel relations are derived
> to obtain the solution.*"

**Atalay 1997, p. 230**:
> "*The singular eigenfunction method provides an analytic treatment
> based on the normal-mode expansion of isotropic scattering,
> monoenergetic transport equation. ... Lately, it has been shown
> that the singular eigenfunction method also serves to evaluate the
> complete spectrum of critical thicknesses and eigenvalues for the
> reflected slabs and spheres with isotropic scattering (Atalay,
> 1995).*"

**Mitsis 1963 (ANL-6787), abstract**:
> "*Transport solutions to the monoenergetic plane, spherical, and
> cylindrical critical problems with isotropic scattering are
> developed by the **method of singular expansion modes**.*"

The two ORPHEUS folders represent **two parametric variations on
the same method**, not two methods:
- `case_method/` (Atalay 1997 family): slab + sphere, **linear
  anisotropy**, half-range bi-orthogonality + first-order Fredholm
  iteration.
- `singular_eigenfunction/` (Westfall-Metcalf 1973 family): cylinder,
  **isotropic scattering**, modified-Bessel-K radial kernels via the
  addition theorem + full-range completeness theorem.

Different geometry, different scattering anisotropy, different
kernel reduction (exponential E_n's vs Bessel K_0's), but
**identical mathematical machinery** above the trusted-library line:
discrete eigenvalue ν₀ (same dispersion function), continuum modes
on (−1, 1), half-range orthogonality, Fredholm reduction of the
boundary condition.

**This is recognised in the ORPHEUS code itself** — both packages
already call `fn_method.core.dispersion.case_nu0` for the dominant
discrete eigenvalue, because the dispersion function is a *medium*
property shared across all Case-method instances regardless of
geometry. See `singular_eigenfunction/__init__.py:39-46`.

### Q2.b — Method-canonical name?

The literature unambiguously uses **"singular eigenfunction
expansion"** (or "singular eigenfunction method" / "singular
expansion modes"):

| Source | Phrase used |
|--------|-------------|
| Case 1960 (foundational) | "elementary solutions" (the *result*); the technique is what later authors named after Case |
| Mitsis 1963 | "method of singular expansion modes" |
| Westfall-Metcalf 1973 (title) | "**Singular Eigenfunction** Solution" |
| Atalay 1997 | "Case's **singular eigenfunction method**" |
| Sanchez 1977 NSE 64 | "**Case singular eigenfunction expansion**" |

**Recommended name**: `singular_eigenfunction/` (KEEP).

`case_method/` violates the user's rule (author-named) AND is not
the standard terminology in the field — "singular eigenfunction
method" is the established phrase.

### Q2.c — Unified folder layout proposal

```
singular_eigenfunction/
├── __init__.py              # rolls up all geometries + anisotropy variants
├── core/                    # shared primitives
│   ├── half_range.py        # was case_method/core/half_range.py
│   ├── x_function.py        # was case_method/core/x_function.py (Atalay X-function)
│   ├── dispersion.py        # was case_method/core/dispersion.py
│   ├── extrapolated_endpoint.py
│   └── bessel_radial.py     # NEW — could host Bessel addition-theorem helpers
│                            #       used by cylinder (currently inline in cylinder/one_group.py)
├── slab/
│   ├── __init__.py
│   ├── one_group.py         # was case_method/slab/one_group.py (Atalay 1997 reflected slab)
│   └── ...
├── sphere/
│   ├── __init__.py
│   ├── one_group.py         # was case_method/sphere/one_group.py (Atalay 1997 reflected sphere via parity flip)
│   └── ...
├── cylinder/
│   ├── __init__.py
│   └── one_group.py         # KEEP (was singular_eigenfunction/cylinder/one_group.py — Westfall-Metcalf 1973)
├── anisotropy/
│   └── linear.py            # was case_method/anisotropy/linear.py (linear-anisotropy half-range relations)
└── origins/
    ├── __init__.py
    ├── slab_sphere.py       # was case_method/origins/derivations.py (Atalay symbolic)
    └── cylinder.py          # was singular_eigenfunction/origins/cylinder_derivations.py (WM symbolic)
```

The `anisotropy/` subfolder is justified because **linear anisotropy
is a transversal feature** — Atalay's slab-sphere parity trick works
for both isotropic and linearly-anisotropic (the half-range
relations are different but the structural reduction to a Fredholm
integral equation is the same). A future `cylinder/linear_anisotropy.py`
could mount on the same `anisotropy/linear.py` half-range machinery
once a benchmark is found.

### Renaming impact note (Q2)

Scope:
- **6 Python files** in `case_method/` to relocate (slab, sphere,
  origins, anisotropy, core).
- **5 test files** in `tests/derivations/test_case_method_*.py`
  (renames: `test_case_method_slab.py` → `test_singular_eigenfunction_slab.py`,
  similarly for sphere/symbolic/parity_flip/z0).
- **1 Sphinx page** `docs/theory/case_method.rst` to merge into
  `docs/theory/singular_eigenfunction.rst` (the singular_eigenfunction
  page already exists and currently covers cylinder only — expanding
  to cover slab+sphere is natural since Atalay explicitly uses the
  parity-flip slab→sphere reduction).
- **1 toctree entry** in `index.rst`.
- **1 import** in `sood_registry/atalay1997.py` (the Atalay case
  registry will keep its filename — atalay1997.py is the *paper*,
  not the *method*; allowed under "author registries OK" rule).
- **1 cross-reference** in `docs/theory/singular_eigenfunction.rst:219`.

The consolidation has **a real architectural payoff** beyond
naming: the X-function machinery (`case_method/core/x_function.py`)
and dispersion machinery (currently shared via `fn_method.core.dispersion`
above the trusted-library line) become first-class members of a
single package. Cross-cutting tests like
`test_carlvik_galerkin_xverif_fn.py` that already cross-check Case-
method results against F_N gain a clearer geometric story.

---

## Q3 — Audit of remaining folders

| Folder | Current | Verdict | Proposed | Rationale | Impact |
|--------|---------|---------|----------|-----------|--------|
| `analytical/` | method-name | KEEP | (same) | "Analytical" is the standard verification-pillar word for closed-form algebraic references (homogeneous-medium k_∞, ratio tests). Could sharpen to `closed_form/` but the existing word is honest. | None |
| `carlvik_galerkin/` | half-author/half-method | RENAME | `galerkin_spectral/` | "Carlvik" is an author. The structural method is **Legendre-Galerkin spectral expansion**. Carlvik introduced *recurrences for matrix elements*; Dahl-Sjostrand 1979 introduced the *Galerkin matrix construction*. Folder should reflect the structural method. The Carlvik recurrences live INSIDE the package (`core/carlvik_recurrences.py` — that filename is OK, "Carlvik recurrences" is a named mathematical object like "Wigner 3-j symbol"). | 5 test files (`test_carlvik_galerkin_*.py`), 1 theory page (`docs/theory/carlvik_galerkin.rst`), index toctree, 6 Python files in package. |
| `case_method/` | author-name | RENAME | merge into `singular_eigenfunction/` | See Q2 above. | Q2 impact note. |
| `cases/` | composition | KEEP | (same) | "Per-method case registries" is honest functional naming. Not method-canonical because it's not a method — it's the case enumeration that composes references. | None |
| `flat_source_cp/` | method-canonical | KEEP | (same) | Standard nomenclature: collision-probability with flat-source-per-region approximation. | None |
| `fn_method/` | method-canonical | KEEP | (same) | F_N method is the established journal-literature name (Siewert et al. 1973–1986). Underscore vs hyphen is a Python-convention choice; the literature writes "F_N" or "F\_N method". | None |
| `mms/` | method-canonical | KEEP | (same) | Method of Manufactured Solutions — universal verification-engineering acronym, not specific to any author. | None |
| `peierls_greens_function/` | author-influenced | RENAME | `trajectory_resolvent/` | See Q1 above. | Q1 impact note. |
| `peierls_nystrom/` | method-canonical | KEEP | (same) | Peierls = the integral form (named after Peierls 1939, but "Peierls integral equation" is now a generic mathematical name like "Boltzmann equation"); Nyström = the discretisation method. Both are textbook method names. **Acceptable — though strictly Peierls is a person, the term is so deeply integrated into transport literature that it functions as a method name.** | None |
| `singular_eigenfunction/` | method-canonical | KEEP & EXPAND | (same) | Will absorb `case_method/`. See Q2. | Q2 impact note. |
| `sood_registry/` | author-name on REGISTRY | KEEP | (same) | User explicitly allowed: "Author names → registry of cases (acceptable)". This is exactly that. | None |

---

## Q4 — Methods not yet implemented (reserve folder names)

Based on local literature (`/workspaces/ORPHEUS/scratch/literature/`)
+ standard transport-theory taxonomy:

| Method | Proposed folder | Purpose | Key reference (LOCAL) | Priority |
|--------|-----------------|---------|------------------------|----------|
| **Spherical-harmonics (P_N)** | `pn_method/` | Modal expansion of angular flux in Legendre polynomials. Multi-region sphere k_eff and flux shapes; the historical alternative to S_N. | Garcia 2017/2019/2021 (3 PDFs LOCAL — sphere shell, exterior of sphere, multi-region sphere). Davison 1957 textbook chapters. | High — Garcia 2021 multi-region sphere is the natural cross-check for trajectory_resolvent multi-region. |
| **Simplified P_N (SP_N)** | `spn_method/` | Asymptotic reduction of P_N in slab geometry, generalised to multi-D. | Makine 2018 (LOCAL — "Exact transport representations of the classical and nonclassical simplified P_N equations"). Brantley-Larsen 2000. | Medium — bridges diffusion and full transport; pedagogically valuable. |
| **Hybrid Galerkin-S_N** | `galerkin_sn_hybrid/` | Morel's collocation-Galerkin angular discretisation as a hybrid between modal and discrete-ordinate. | Morel 1989 (LOCAL — "A Hybrid Collocation-Galerkin-Sn Method for Solving the Boltzmann Transport Equation"). | Low — research-grade, narrow application. |
| **B_N method** | `bn_method/` | Sister of F_N: collocates F_N on the *boundary* points instead of interior; used for diffuse-reflection benchmarks. | Sood/Forster/Parsons 2003 cite it; Garibba-Rojas 1980s (NOT local). | Low — niche; F_N covers most use-cases. |
| **Method of characteristics (MoC)** | (already in `discrete/`) | NOT a continuous reference — discrete by construction. | Askew 1972 (LOCAL — "A Characteristics Formulation of the Neutron Transport Equation"), Halsall 1980 (LOCAL — CACTUS), Mazumdar 2015 (LOCAL — review). | N/A — already covered as a discrete method. |
| **Integral transport (Pomraning-Siewert family)** | (subsumed by `trajectory_resolvent/` after Q1 rename) | The pure integral-equation form for sphere/slab without trajectory tracking; PS-1982 and Sanchez 1986 are the canonical specular-BC references. | Pomraning-Siewert 1982 (LOCAL), Sanchez 1986 (LOCAL), Pomraning 1989 (LOCAL — "The Transport Equation in General Geometry"). | Already covered — PS-1982 verification is in `tests/derivations/test_peierls_greens_function_xverif_ps1982.py`. |
| **Escape-probability method** | `escape_probability/` | Wigner-Seitz cell probabilities; the resonance-self-shielding ancestor of CP. | Carlvik 1965 Geneva (cited by Sanchez 1977 — NOT local), Carlvik 1967 NSE 30 (LOCAL — finite cylinder/cuboid CP via E_3 kernels). Benoist 1981 (LOCAL — "Integral Transport Theory Formalism for Diffusion Coefficient Calculations in Wigner-Seitz Cell"). | Medium — used heavily in lattice physics; related to but distinct from `flat_source_cp/`. |
| **Reed's method** (asymptotic / boundary-layer corrections) | (no clear single-name candidate) | Reed-Hill thermal-Maxwell-Boltzmann correction (1971, NSE 46) is sometimes called "Reed's method" — but it's a *boundary correction*, not a discretization method. Doesn't warrant a folder. | Reed-Hill 1971 NSE — NOT local. | None — not a structural method. |
| **Method of images** (slab) | `method_of_images/` *(maybe)* | Eigenfunction expansion via reflection symmetry — used as cross-check for slab problems with reflective BC. Pedagogically valuable but rarely a primary verification method. | Case-Zweifel 1967 ch.3, Hetrick 1971 textbook. | Low — most modern slab benchmarks use F_N or singular_eigenfunction directly. |
| **Spectral collocation (Chebyshev / Gegenbauer)** | `spectral_collocation/` | Chebyshev / Gegenbauer-collocation methods for transport — recent work on stability for high-anisotropy. | Garcia-Siewert 1979 / 1980s (NOT all local). Hauser-Pomraning 1980s. | Low — niche; most useful as Branch-1 (symbolic verification). |

**Strong recommendation**: create empty `pn_method/` folder NOW with
`__init__.py` stub citing Garcia 2021. Garcia 2021 is the multi-region
sphere P_N that ORPHEUS already cross-checks against in
`test_peierls_greens_function_garcia2021.py`. Having `pn_method/`
exist as a reserved name (even empty) makes the cross-check's
target-vs-reference asymmetry visible structurally — currently
ORPHEUS uses Garcia 2021 as a *truth set without a method-of-record*,
which is a hidden coupling.

---

## Literature consulted

All sources LOCAL in `/workspaces/ORPHEUS/scratch/literature/`:

1. **Sanchez, R., Mao, L., Santandrea, S. (2002)**. "Treatment of
   Boundary Conditions in Trajectory-Based Deterministic Transport
   Methods." *Nuclear Science and Engineering* **140**, 23–50.
   [Read pages 1-10; Eq. (15) periodic-trajectory closure verified.]
2. **Sanchez, R. (1986)**. "Integral form of the equation of transfer
   for a homogeneous sphere with linearly anisotropic scattering."
   *Transport Theory and Statistical Physics* **15**(3), 333–343.
   DOI: 10.1080/00411458608210456.
   [Full Appendix read; Eq. (A4) is THE resolvent ORPHEUS implements.]
3. **Pomraning, G.C., Siewert, C.E. (1982)**. "On the integral form
   of the equation of transfer for a homogeneous sphere." *J.
   Quant. Spectrosc. Radiat. Transfer* **28**(6), 503–506.
   [Full 4-page note read; Eq. (14) `T(µ) = [1 − α e^{−2Rµ}]^{−1}`
   is the structurally-independent precursor.]
4. **Westfall, R.M., Metcalf, D.R. (1973)**. "Singular Eigenfunction
   Solution of the Monoenergetic Neutron Transport Equation for
   Finite Radially Reflected Critical Cylinders." *Nuclear Science
   and Engineering* **52**, 1–11. DOI: 10.13182/NSE73-A23285.
   [Full 11 pp + Appendix read; explicit naming "singular eigenfunction
   expansion technique by Case", confirming Q2 method-class identity.]
5. **Atalay, M.A. (1997)**. "The reflected slab and sphere criticality
   problem with anisotropic scattering in one-speed neutron transport
   theory." *Progress in Nuclear Energy* **31**(3), 229–252.
   DOI: 10.1016/0149-1970(95)00094-1.
   [First 10 pp read; abstract + p. 230 explicit "Case's singular
   eigenfunction method".]
6. **Mitsis, G.J. (1963)**. "Transport Solutions to the Monoenergetic
   Critical Problems." Argonne National Laboratory report ANL-6787
   (PhD thesis, U. Michigan).
   [TOC + abstract + first chapter read; "method of singular
   expansion modes" verbatim — pre-dates Westfall-Metcalf and is the
   foundational source they both cite.]

Additional local sources surveyed (not deeply read for this memo,
but contents confirmed via package docstrings + abstracts):

7. Sanchez 1977 NSE 64 (CP collision probability) — Class-B IC paper.
8. Sanchez 1982 (NSE 80) — review of neutron transport approximations.
9. Carlvik 1967 NSE 30 — finite cylinder/cuboid CP via E_3 kernels.
10. Garcia 2017/2019/2021 — P_N family for sphere geometries.
11. Makine 2018 — exact transport SP_N representations.
12. Morel 1989 — Galerkin-S_N hybrid.
13. Ligou 1982 ch.8, Stammler 1983 ch.4+6, Stacey 2007 ch.9, Hebert
    2009 ch.3 — textbook chapters (all LOCAL).
14. Sood/Forster/Parsons 1999/2003 — LA-13511 benchmark catalogue.
15. Kaper-Lindeman-Leaf 1974 — slab+sphere F_N benchmark values.
16. Siewert-Thomas 1986 — 2G F_N.
17. Pomraning 1989 — Transport Equation in General Geometry.

---

## Literature requested from user

**None.** Every claim in this memo is supported by a paper already in
`/workspaces/ORPHEUS/scratch/literature/`. The user's local archive
covers the entire chain Mitsis 1963 → Case 1960 (cited but not
needed locally for the renaming claim) → Westfall-Metcalf 1973 →
Sanchez 1986 → Sanchez 2002 → Atalay 1997, and the parallel
Pomraning-Siewert 1982 + Sanchez 1977 + LA-13511 verification
chain.

If the user wants to validate any claim deeper than this memo
asserts, the four candidates worth reading PDF-fully would be:

1. **Case 1960 *Ann. Phys. NY* 9, 1** — to anchor the "Case method"
   foundational claim. NOT local. Worth requesting if you want a
   primary-source citation rather than secondary citations from
   Mitsis/Westfall-Metcalf/Atalay.
2. **Carlvik 1968 NSE 31, 295** — the original Carlvik recurrences.
   NOT local; only the 1967 NSE 30 finite-cuboid CP paper is local.
   Would harden the Q3 `carlvik_galerkin/` rename rationale.
3. **Dahl-Sjostrand 1979** — already cited extensively in
   `carlvik_galerkin/__init__.py` but the 1979 NSE 69 paper itself
   is NOT local. The Dahl-Sjostrand 1979 anisotropic eigenvalue paper
   IS local — that may be sufficient depending on the user's intent.
4. **Smith 1969 PhD thesis** (NCSU) — Westfall-Metcalf cites this
   for the periodic-trajectory closure that pre-dates everything.
   Likely paywalled/scarce; only worth chasing if the historical
   priority on `trajectory_resolvent/` matters.

None of these block the renaming work. Recommendation: proceed with
the renames using the evidence in this memo; if Case 1960 turns out
to be needed for a Sphinx narrative section, add it then.
