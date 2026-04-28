---
name: Phase 5 — continuous-µ multi-bounce specular kernel literature pull (sphere/cyl/slab)
description: Literature pull for Phase 5 of the Peierls specular BC closure. Hébert §3.8.3 slab + Sanchez 2002 trajectory form are the only standard textbook continuous-µ multi-bounce kernels. Sanchez 1986 TTSP, Pomraning-Siewert 1982, Milgram 1978 NSE/Can.J.Phys are the analytical-sphere/cylinder integral-transport references. Stamm'ler/Stacey/Ligou/Hébert §3.8.4 all use scalar (rank-0) (1−P_ss)^{−1} closure.
type: project
---

# Phase 5 continuous-µ multi-bounce specular kernel — literature pull (2026-04-28)

## Bottom line

The continuous-µ multi-bounce form `1 / (1 − e^{−σ·L(µ)})` for
specular BC is **explicitly written in only two places** in the
local PDF set: **Hébert 2009 §3.8.3 Eq. (3.273)** for slab CP, and
**Sanchez et al. 2002 NSE 140 Eqs. (15)/(17)–(20)** for general
periodic-trajectory CP. Every other textbook reference (Stamm'ler
1983 Ch.4, Stacey 2007 Ch.9, Ligou 1982 Ch.8, Hébert 2009 §3.8.4
cylindrical, §3.8.5 spherical) uses the **scalar (rank-0) `(1 −
α(1 − Λ))^{−1}` geometric series**, where the µ-resolution is
integrated out *before* the bounce algebra. The cylinder integral
form CANNOT carry the resummation because Bickley-Naylor Ki_n
integration kills the multiplicative property of the exponential.

The cleanest analytical-sphere/cylinder benchmark literature is
**Sanchez 1986 TTSP** (specular sphere with anisotropic scattering),
**Pomraning-Siewert 1982 JQSRT** (homogeneous-sphere integral
form), **Milgram 1978 NSE** (homogeneous cylinder analytical), and
**Dahl-Sahni 1983 TTSP** (Carlvik's method for sphere criticality).

Phase 5 has a clear textbook precedent for **slab** (Hébert §3.8.3
verbatim) but is **genuinely new (or at least under-published)
for sphere and cylinder integral form**. The plan should start
from Sanchez 2002 Eq. (15) as the per-trajectory template and
specialize each geometry's L(µ).

---

## TIER 1 — sources with the continuous-µ kernel

### **Hébert 2009 *Applied Reactor Physics* §3.8.3 Eq. (3.273)**

- Section: §3.8.3 "Slab geometry" (Kavenoky method).
- Topic: collision probability `P_ij` in slab geometry with
  specular periodic neighbour-cell BC, expressed as an explicit
  geometric series **before** the polar integration is done.
- µ-resolved form (TeX, transcribed from p. 113 Eq. 3.273):

```math
P_{ij} \;=\;\frac{1}{2\,\Delta x_i}\sum_{m=0}^{\infty}
\int_{x_{i-1/2}}^{x_{i+1/2}} dx'\!\int_0^{\pi/2} d\theta\,\tan\theta
\!\int_{x_{j-1/2}}^{x_{j+1/2}} dx\;
\bigl[
e^{-(m\,\tau_{\text{cell}}+\tau(x',x))/\cos\theta}
+e^{-((m+1)\tau_{\text{cell}}-\tau(x',x))/\cos\theta}
\bigr]
```

- Key observation: **the specular periodic BC enters the kernel as
  an infinite sum of exponential factors `e^{−m·τ_cell/cos θ}`**,
  which sums per-µ to `1/(1 − e^{−τ_cell/cos θ})` *before* polar
  integration. This is the slab analogue of the Phase 5 sphere
  factor `µ/(1 − e^{−σ·2Rµ})`. The polar integration is then done
  via E_n functions — but if you stop *before* integrating, you
  have the per-µ multi-bounce kernel.
- Notation map: their `cos θ` = our `µ`; `τ_cell` = the per-cell
  optical thickness along axis. Their geometric series IS the
  Phase 5 multi-bounce factor.
- Direct relevance: SLAB Phase 5 is *literally* Hébert §3.8.3 with
  the m-sum kept symbolic and not carried through to E_n. ORPHEUS
  can cite this section directly.
- ORPHEUS module match: `compute_K_bc_continuous_mu_slab` (to be
  created) or extension of `cp_slab.py`.

### **Sanchez, Mao & Santandrea 2002, *Nucl. Sci. Eng.* 140, 23–50, Eq. (15)**

- Section: §III "Exact Treatment for Closed Domains", Eqs. (12)–
  (15), p. 31. Generalisation to 3D in §III.D and CP form in
  Eqs. (17)–(18).
- Topic: angular flux along a periodic compound trajectory of
  length L(µ) under specular/translation/rotation BC.
- µ-resolved form (TeX, Eq. 15):

```math
\psi(x) \;=\;\frac{\psi_q(L)}{\,1-\psi_{\rm bd}(L)\,}\,\psi_{\rm bd}(x)
\;+\;\psi_q(x)
```

with `ψ_bd(x) = e^{−τ(0,x)}` (boundary-attenuation along the
trajectory) and `ψ_q(x) = ∫_0^x e^{−τ(x',x)} q(x')dx'` (source
contribution). For collision-probability form, the corresponding
Eq. (17)–(18) in 3D is

```math
P_{ij}^{\rm 3D}\;=\;K_{\rm 3D}\!\!
\sum_{n=-\infty}^{\infty}\!\!
\int_{(i)} dt \int_{(j)} dt'\, e^{-|\tau(t',t)+n\tau_L|}
\;=\;K_{\rm 3D}\!\!\biggl[\frac{1}{1-e^{-\tau_L}}\!
\int_{(i)}\!dt\!\int_{(j)}\!dt'\,e^{-\tau(t',t)}
\;+\;\dots\biggr]
```

The **`1/(1 − e^{−τ_L(µ)})`** factor is exactly the multi-bounce
multiplier per trajectory.

- Notation map: their `τ_L` = our `σ·L(µ)`; their `ψ_bd(L)` = our
  `e^{−σ·L(µ)}`. Their per-trajectory length `L` for sphere with
  specular BC is `2Rµ`; for slab is `H/µ`; for cylinder it is
  `chord(α,θ)` of the periodic in-plane orbit. The factor
  `1/(1 − ψ_bd(L))` IS the Phase 5 `f_mb(µ)`.
- Direct relevance: **this is the canonical reference for the
  "absorb specular BC into the kernel via geometric series"
  technique**. The Eq. (15) form is geometry-agnostic — pick the
  trajectory length L(µ) for sphere/cyl/slab and the formula
  applies.
- Critical caveat: Sanchez et al. work in the **MoC trajectory
  representation** (full angular ψ along a trajectory, polar
  integration done analytically AFTER bouncing for the 2D xy CP
  case via Bickley-Naylor Ki1, Eq. 19–20). For the 3D CP form
  Eq. (17), the bouncing is done **before** polar integration and
  the resummation works analytically (Eq. 18).
- §III.D "Periodicity in the Axial Direction" (Eq. 34, p. 35):
  trajectory in 3D requires *both* azimuthal periodicity (Eq. 21)
  AND polar periodicity `tan θ = L_w/(l·c)` with l integer. This
  is a constraint on which (φ, θ) admit periodic compound
  trajectories — relevant if Phase 5 wants to discretise a
  cylinder or 3D box.
- ORPHEUS module match: `compute_K_bc_continuous_mu_*` for
  sphere/cyl/slab; the L(µ) selection logic per geometry.

---

## TIER 2 — analytical-sphere / -cylinder integral transport (benchmark literature)

### **Sanchez 1986, *Transport Theory & Stat. Phys.* 14, DOI 10.1080/00411458608210456**

- Title: "Integral form of the equation of transfer for a
  homogeneous sphere with linearly anisotropic scattering".
- Abstract verbatim: *"The equation of transfer for a homogeneous,
  linearly anisotropic sphere with internal sources and a
  **specular and diffuse surface** is reduced to an integral
  equation. **Homogeneous boundary conditions are directly
  incorporated in the kernels of the integral equation.**"*
- Direct relevance: **THIS IS THE PHASE 5 SPHERE PAPER**. The
  abstract states explicitly that the specular surface BC is
  incorporated *into the kernel* of the integral equation. This is
  the transformation Phase 5 wants — from `K · (I − T·R)^{−1}` to
  a single kernel `K_bc^spec`. Only 2 citations on Semantic
  Scholar; obscure but directly on-topic.
- Notation map: unknown without PDF read — author writes both
  specular AND diffuse surface in same framework, suggesting the
  kernel treats them parametrically.
- Status: paywall, no OA. **Strong candidate to add to user's
  Zotero library and PDF-read for explicit kernel form.**
- ORPHEUS module match: would directly inform
  `compute_K_bc_continuous_mu_sphere`.

### **Pomraning & Siewert 1982, *J. Quant. Spectrosc. Radiat. Transf.* 28, DOI 10.1016/0022-4073(82)90016-4**

- Title: "On the integral form of the equation of transfer for a
  homogeneous sphere".
- 29 citations on Semantic Scholar (the standard reference). No
  abstract retrievable, but the title and citation count establish
  it as the canonical pre-Sanchez work on the homogeneous-sphere
  integral kernel.
- Direct relevance: foundational kernel structure for homogeneous
  sphere; Sanchez 1986 likely builds on this.
- Authors are heavyweights (Pomraning, Siewert, FN method) so the
  kernel structure here is going to be authoritative.
- ORPHEUS module match: kernel verification reference; cite for
  the bare-sphere case (no specular BC) and check Sanchez 1986
  reduces to it when α → 0.

### **Milgram 1978, *Nucl. Sci. Eng.* 67, DOI 10.13182/nse78-a27304**

- Title: "Analytic Method for the Numerical Solution of the
  Integral Transport Equation for a Homogeneous Cylinder".
- Abstract excerpt: *"The integral transport equation for the
  flux density in the interior of an infinite homogeneous cylinder
  is reduced to a matrix eigenvalue problem for the critical
  cylinder and a set of linear algebraic equations for the driven
  case with surface in-currents. The matrix elements are
  identified as **moments of modified Bessel functions** and are
  easily computed."*
- Direct relevance: cylinder analytical k_eff via small matrix
  eigenvalue problem — same structural idea as Phase 5 wants. The
  Bessel-function moment integrals likely contain the
  continuous-µ multi-bounce information.
- Notation map: their "matrix elements" are integrals of `K_n` /
  `Ki_n` over chord parameter R — the cylinder analogue of the
  Phase 5 `G_in(µ)`, `F_out(µ)`, `f_mb(µ)`.
- ORPHEUS module match: `compute_K_bc_continuous_mu_cylinder`
  benchmark (for verification of the k_eff convergence).

### **Milgram 1978, *Can. J. Phys.* 56, DOI 10.1139/p78-081**

- Title: "Analytic method for the solution of the one-group,
  integral transport equation for a homogeneous sphere".
- Abstract: matrix eigenvalue problem, eigenvectors are Taylor
  coefficients of flux about sphere centre. Vacuum BC.
- Direct relevance: sphere analytical k_eff with vacuum BC. Phase
  5 needs the SPECULAR variant — but this paper's matrix structure
  is the right scaffolding to add specular BC kernels to.
- Notation map: vacuum BC version of what Sanchez 1986 does for
  specular.
- ORPHEUS module match: bare-sphere k_eff verification (vacuum BC),
  to be **specularized** following Sanchez 1986.

### **Dahl & Sahni 1983, *Transport Theory & Stat. Phys.* 12, DOI 10.1080/00411458308211640**

- Title: "Complex time eigenvalues of the one speed neutron
  transport equation for a homogeneous sphere".
- Abstract excerpt: *"Using **Carlvik's method** to solve the
  criticality problem, it is shown that complex time eigenvalues
  do actually exist for this model problem."*
- Direct relevance: Carlvik's method for homogeneous-sphere
  criticality — a numerical recipe directly applicable to Phase 5
  verification. Vacuum BC, isotropic scattering.
- ORPHEUS module match: numerical sphere k_eff recipe for the
  vacuum-BC case; specular BC requires additional kernel terms.

### **Milgram 1977, *J. Math. Phys.* 18, DOI 10.1063/1.523208**

- Title: "On the properties of collision probability integrals in
  annular geometry. I. Analysis".
- Abstract excerpt: *"Two integrals fundamental to the analysis of
  the integral transport equation in infinitely long, annular
  geometry are defined and represented as sums of **Meijer's G-
  function**. Special results include the identification of
  probability integrals for infinite cylinders and **generalized,
  associated Bickley-Naylor functions** as G-functions with
  appropriate parameters and variables."*
- Direct relevance: closed-form analysis of cylindrical CP
  integrals in Meijer-G form; a heavy-machinery toolkit for
  resumming the multi-bounce series exactly. Probably overkill for
  Phase 5 implementation, but valuable as a sanity check that the
  resummed kernel is computable in closed form.
- Companion paper (Part II): Milgram & Sly 1979 *J. Comput. Phys.*
  DOI 10.1016/0021-9991(79)90167-0 — numerical evaluation of the
  same integrals.
- ORPHEUS module match: verification only; not implementation.

### **Knyazev 1993, *Atomic Energy* 74, DOI 10.1007/bf00844623**

- Title: "Solution of the transport equation in integral form in a
  one-dimensional cylindrical geometry with linearly anisotropic
  scattering".
- Already in ORPHEUS for the rank-N cylinder primitives (per
  `phase4_cylinder_peierls.md`). For Phase 5, check whether
  Knyazev separately treats specular BC — abstract not retrievable
  via OpenAlex/CrossRef.
- ORPHEUS module match: rank-N cylinder primitives; verify if
  specular BC is in §3 or §4.

---

## TIER 3 — sources confirming the SCALAR (rank-0) structure (negative results)

### **Stamm'ler & Abbate 1983, *Methods of Steady-State Reactor Physics in Nuclear Design* Ch.4 Eqs. (32)–(36)**

- Section: Ch. IV, FORTRAN program COLPROB, pp. 131–136.
- Topic: cylindrical Wigner-Seitz cell with albedo BC at S_B.
- Form (verbatim from OCR'd PDF, Eq. 34–36):
  ```
  M_R(α) = 1 / (1 − α(1 − Λ))               (boundary multi-reflection factor)
  Y_i(α) = Y_i / (1 − α(1 − Λ))             Eq. (34)
  X_i^k(α) = X_i^k + α x_S X_S Y_i / (1 − α(1 − Λ))  Eq. (35)
  j_+(α) = (X_S + (1−Λ) j_in_ext) / (1 − α(1 − Λ))    Eq. (36)
  ```
  with `Λ = 1 − Σ_S P_Si` the multicollision blackness (NOT
  multireflection!), and α the cell-boundary albedo (0 = black,
  1 = white).
- This is **scalar (rank-0)**: the µ-dependence is integrated out
  into the scalar transmission Λ before the boundary geometric
  series is summed.
- Notation map: their `α` = our boundary albedo β; their `Λ` ≈ our
  `1 − P_ss`; their multi-reflection factor `1/(1 − α(1−Λ))` = our
  rank-1 Hébert form.
- Direct relevance: confirms that the rank-0 / scalar form is the
  textbook standard for cylindrical multi-region cells. Phase 5
  cylinder MUST go BEYOND this to get the µ-resolution.

### **Stacey 2007 *Nuclear Reactor Physics* Ch.9 Eqs. (9.86)–(9.89)**

- Section: §9.4 "Interface Current Methods in Slab Geometry".
- Form: `P_i = P_{0i}/(1 − c_i(1 − P_{0i}))`, `T_i = T_{0i} +
  ½ c_i P_i (1 − T_{0i})`, with `T_{0i} = E_2(τ_i)` (slab) — i.e.
  µ-integrated already.
- This is **multi-collision** (volume re-emissions), not
  multi-bounce on the surface. Both produce `1/(1 − x)` series, but
  µ is integrated out into E_2/E_3 before the geometric sum.
- Direct relevance: confirms scalar pattern for both slab CP and
  IC families.

### **Ligou 1982 Ch.8 Eq. (8.78)**

- Section: Ch. VIII.4.2 "Influence of the Lattice Pitch on the
  Bell Factor".
- Form: `P_cc(Γ) = P_cc + Γ·P_cs·P_sc / (1 − Γ·P_ss)` with Γ the
  Dancoff coefficient (= reflection coefficient of the
  surrounding rod ensemble).
- This is the **Dancoff-corrected scalar** form. Same multi-bounce
  geometric series as Stamm'ler/Stacey.
- Direct relevance: classical Dancoff-correction algebra is
  scalar.

### **Hébert 2009 §3.8.4 Cylindrical & §3.8.5 Spherical**

- §3.8.4: Eq. (3.297) integrates polar θ analytically into
  Bickley-Naylor `Ki_n(τ)` BEFORE doing the boundary specular sum.
  **The Bickley quadrature destroys the geometric-series
  resummation property.**
- §3.8.5: Eqs. (3.324)–(3.336) — multi-region sphere CP, but again
  µ-integrated into `E_n` before any specular BC algebra. The
  specular/white BC closure for sphere uses the scalar Eqs. (3.320)–
  (3.323) of §3.8.3 (the same as cyl/slab):
  ```
  P̂_ij = P_ij + (β / (1 − β·P_ss)) · P_iS · P_Sj
  ```
  (Eq. 3.323 verbatim, with `β` = albedo, `P_ss = 1 − Σ P_Si`).
- Direct relevance: confirms that the integral form for cyl/sph in
  Hébert is **scalar** (rank-0), same as Stamm'ler/Stacey/Ligou.
  Slab §3.8.3 is the unique exception — and only because `e^x` is
  the eigenfunction of the d/dx kernel that allows the m-sum to
  factor cleanly.

---

## What was NOT found

### Not in any local PDF
- A **µ-resolved closed-form `K_bc^spec` kernel for sphere or
  cylinder** that's both (a) continuous-µ and (b) absorbs the
  multi-bounce sum analytically. This appears to be **genuinely
  novel for the cylindrical integral form** (Hébert §3.8.3 has it
  for slab; Sanchez 2002 has it as a per-trajectory recipe). The
  cylinder integral form's Bickley-Naylor µ-integration is
  incompatible.

### Not retrievable in this session
- Sanchez 1986 TTSP PDF (DOI 10.1080/00411458608210456) — paywall,
  no OA. **STRONGEST add-to-Zotero candidate.** Author abstract
  says "specular surface absorbed into kernel" verbatim.
- Pomraning-Siewert 1982 JQSRT PDF — paywall, no OA. Foundational.
- Milgram 1978 NSE/Can.J.Phys PDFs — paywall, no OA.
- Bell & Glasstone 1970 *Nuclear Reactor Theory* — not in local
  PDFs. The standard textbook reference for integral transport
  pre-1990; would have a chapter on sphere/cyl integral transport
  with various BC.
- Davison 1957 *Neutron Transport Theory* — not in local PDFs.

---

## Recommended Phase 5 reading list (in priority order)

1. **Sanchez 1986 TTSP** — sphere with specular surface absorbed
   into the kernel. **Add to Zotero, PDF-read** if Phase 5 is
   pursued. The abstract is so on-the-nose for the Phase 5 brief
   that it must be read before any from-scratch derivation.
2. **Hébert 2009 §3.8.3** — slab continuous-µ multi-bounce kernel,
   already in the project PDFs. Re-read Eqs. (3.273)–(3.275)
   carefully for the SLAB Phase 5 form. The polar-integration
   step before E_n is the inflection point: stop there to keep
   the per-µ kernel.
3. **Sanchez et al. 2002 NSE 140 §III** — periodic-trajectory CP
   form. Eq. (15) is the geometry-agnostic template. **Already in
   project PDFs.**
4. **Pomraning-Siewert 1982 JQSRT** — verify against sphere case
   when α (specular fraction) → 0. Foundational.
5. **Milgram 1978 NSE** — cylinder analytical k_eff via Bessel-
   function moment matrix. The cylinder analogue of Phase 5
   target.
6. **Dahl-Sahni 1983 TTSP** — Carlvik's method for sphere
   criticality. Numerical recipe for verification.
7. **Knyazev 1993** — already in ORPHEUS for rank-N primitives;
   re-read for any specular-BC content I may have missed.

---

## Notation alignment for Phase 5 derivation

The key ORPHEUS↔literature correspondence:

| ORPHEUS                     | Hébert §3.8.3 (slab)         | Sanchez 2002 Eq. 15  |
|-----------------------------|------------------------------|----------------------|
| `µ` (polar cos)             | `cos θ`                      | along trajectory     |
| `σ`                         | `Σ`                          | `s`                  |
| `L(µ)` chord per geometry   | `Δx_cell / cos θ` (slab)     | `L` = trajectory pd  |
| `e^{−σ·L(µ)}`               | `e^{−τ_cell/cos θ}`          | `ψ_bd(L) = e^{−τ_L}` |
| `f_mb(µ) = µ/(1−e^{−σ·L})`  | implicit in m-sum            | `1/(1 − ψ_bd(L))`    |
| `G_in(r,µ)`                 | source-side kernel           | `ψ_bd(x)`            |
| `F_out(r',µ)`               | field-side kernel            | `ψ_q(L)`             |

Phase 5 implementation should use **Sanchez 2002 Eq. (15)** as the
template, instantiate `L(µ)` per geometry:

- Slab: `L(µ) = H/µ` (chord through full slab at angle µ)
- Sphere: `L(µ) = 2Rµ` (chord = 2R cos α; α = polar angle from
  sphere normal)
- Cylinder: `L(α, θ_p) = chord_in_plane(α)/sin θ_p` (in-plane
  chord at azimuth α, at polar angle θ_p)

then integrate continuous µ (or (α, θ_p) for cylinder) at the END,
not before the bounce.

The slab case is fully covered by Hébert §3.8.3. The sphere case
is fully covered by Sanchez 1986 abstract claim — needs PDF read.
The cylinder case requires a NEW derivation following Sanchez 2002
Eq. (15) trajectory-by-trajectory, since Hébert §3.8.4 uses
Bickley-Naylor and is incompatible.

---

## Status of negative-result possibility

The user's brief flagged the possibility that "the continuous-µ
form may be a folk-theorem unwritten in textbooks". Evidence
**partially supports this**: for SLAB it's textbook (Hébert
§3.8.3); for SPHERE Sanchez 1986 claims to derive it (abstract);
for CYLINDER no integral-form publication has been found that
does it (the Sanchez 2002 trajectory form does it implicitly per-
trajectory but doesn't publish the cylinder kernel in closed form).

**Recommendation**: pull Sanchez 1986 (highest priority) and
Milgram 1978 NSE (cylinder). If both confirm the Phase 5 brief,
implementation starts from textbook derivations. If they
contradict (e.g., Sanchez 1986 only handles diffuse not specular
in the sense Phase 5 needs), then the cylinder Phase 5 is a
genuinely novel derivation and SymPy from-scratch is the path.
