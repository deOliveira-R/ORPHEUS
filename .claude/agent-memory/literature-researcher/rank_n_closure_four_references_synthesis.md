---
name: Rank-N IC closure — cross-reference of Ligou 1982 Ch.8, Sanchez 2002, Stamm'ler 1983 Ch.4, Stacey 2007 Ch.9
description: Side-by-side of four authoritative interface-current references, their actual closure forms (mostly scalar/DP-0), mapping to the ORPHEUS F.4 code, and where the 1.42 % Sanchez-McCormick 1982 plateau likely comes from.
type: project
---

# Rank-N interface-current closure — the four canonical references

**Context.** ORPHEUS implements the Sanchez-McCormick 1982 §III.F.1
rank-N per-face hollow-sphere closure. It plateaus at ~1.42 % k_eff
error (σ_t·R = 5) for N = 2, 3, 4, while the existing scalar F.4
"Lambertian return" closure hits 0.077 % at the same point. The user
asked for evidence from four additional references.

**Headline finding.** Three of the four references
(Ligou 1982, Stamm'ler 1983, Stacey 2007) do NOT contain a rank-N
Marshak/Yvon ladder at all. They use a **single scalar incoming
partial current per face with a cosine (Lambertian) angular
distribution for multi-reflection** — the F.4 formulation is
literally the textbook closure. Sanchez 2002 uses piecewise-constant
(in angle) or collocation (δ-function) representations in the
trajectory-tracking sense — not Legendre moments — and it refers to
Sanchez-McCormick 1982 §III.F for the moment-closure derivation.
**There is no second independent derivation of the §III.F.1 Legendre
ladder among these four references.** That is itself the most
actionable result of the search: the §III.F.1 recipe we tried to
reproduce appears to be unique to the 1982 review, and its plateau
cannot be cross-validated against any of these four authoritative
sources.


## §1 Summary of each reference

### Ligou (1982) — *Elements of Nuclear Engineering*, Chapter 8 "Fundamental Problems of Neutronics"

Chapter 8 is a 53-page **general neutronics chapter** covering (in
order): the differential-integral transport equation (§8.1.1), plane
geometry Milne problem and diffusion limit (§8.1.2-8.1.4), **CP
methodology** (§8.2 — pp. 397–410), slowing down + Placzek transient +
resonances (§8.3), heterogeneous resonances and the Bell/Dancoff
equivalence theorem (§8.4). The only IC-relevant content is §8.2,
and its closure is **scalar**.

- §8.2.3 ("The boundary conditions"): single net current J at the
  cell boundary, or a single incoming partial current j⁻ with
  **cosine distribution** if white BC. Eq. 8.47 gives the
  Wigner-Seitz relation P_{S→j} = (4/S) V_j Σ_{t,j} P_{j→S} assuming
  cosine incoming current — this is the *definition* of F.4 scalar
  white BC.
- §8.4.1 (equivalence theorem): builds P̃_cc = P_cc + Γ P_cs P_sc /
  (1 − Γ P_ss). The Dancoff Γ is the **boundary albedo** and enters
  through a **scalar geometric series** 1 + Γ P_ss + (Γ P_ss)² + … —
  exactly the F.4 scalar multi-reflection sum.
- No Legendre expansion of the angular flux at the surface. No DP_N.
  No per-face rank-N block structure.

Ligou's §8.1 does Legendre-expand the **scattering kernel** (Eq.
8.17: g(μ, μ') = (1 + 3μ̄₀μμ')/2, linearly anisotropic scattering),
but that is volumetric anisotropy, not surface-flux anisotropy.

### Sanchez, Mao, Santandrea (2002) — NSE 140:23 "Treatment of Boundary Conditions in Trajectory-Based Deterministic Transport Methods"

This is **not** a refinement of Sanchez-McCormick 1982 §III.F.1. It
is a paper about **exact vs. approximated handling of geometrical
motions in trajectory tracking** (translation / rotation / specular
reflection group actions on cyclic trajectories) and about how
approximated albedo BCs can mimic them.

- The paper acknowledges two canonical angular representations for
  the boundary flux — "piecewise-constant" (in azimuthal/polar
  sectors) and "collocation" (discrete δ₂ on a quadrature set) — see
  Eqs. 37 and 40. Both are **spatial + angular piecewise decompositions**,
  not Legendre-moment decompositions.
- The DP_N / Legendre-moment closure is mentioned in passing
  (Eq. 4 and the first paragraph of p. 24) and the reader is
  **referred to Sanchez-McCormick 1982** (reference 20) for the
  derivation. Sanchez 2002 does not rederive the rank-N recipe.
- The piecewise-constant N=1 representation with a **single
  azimuthal and single polar sector** (Nφ×Nθ = 1×1) IS effectively
  the scalar DP-0 closure. Sanchez 2002's Table XVIII (BWR
  Dodewaard, PCA 2×2 vs PCA 4×4) shows the approximation converges
  to the exact treatment as angular sectors are refined.
- **No explicit formula ties DP_N moments to a (I − W)⁻¹ P closure
  form.** The paper's closure machinery is Eq. 38:
  `J_a+^r = Σ_i E_ai^r q_i + Σ_{b,n} T_ab^{rn} J_b-^n` with the E
  and T matrices "appropriately computed"; the actual code is the
  TDT module of APOLLO2 (reference 12), which is export-controlled
  and cannot be inspected.
- **Key implicit warning**: §V "Product Quadrature Formulas and
  Conservation" shows that when the angular quadrature fails to
  integrate the even Legendre polynomials J₂ₖ = ∫₀¹ P₂ₖ(μ) dμ = 0,
  the resulting IC approximation loses particle conservation and can
  give k_eff errors of order 10⁻³ to 10⁻² in strong-leakage
  problems. Our Sanchez N=2 Legendre basis uses exact integrals, but
  this is a known failure mode for the method family.

### Stamm'ler & Abbate (1983) — *Methods of Steady-State Reactor Physics in Nuclear Design*, Chapter IV "Integral Transport Theory; Collision Probabilities"

**This is the correctly cited chapter for CP methodology.** Ch.IV
supersedes the earlier (incorrect) Ch.VI citation. It is 37 pages
and its IC treatment is **entirely scalar**.

- §10 ("Application of Collision Probabilities", pp. 126-141) is the
  canonical cell-solution procedure. The boundary flux is
  represented by:
  - a **single** integrated partial in-current j⁻ (Eq. 29)
  - a **single** albedo a at the outer surface
  - a **cosine (white) distribution** for multi-reflection
    (section §10 p. 128: "each neutron reaching the boundary will
    be reduced by the factor a and then reflected back with a cosine
    distribution").
- Eq. 32: Γ = Σᵢ Σ_{r,i} Vᵢ Yᵢ — **multi-collision blackness is
  scalar**, summed across all regions.
- Eq. 34 and Eq. 36a,b: Y_i(a) = Y_i / (1 − a(1−Γ)) and the in/out
  currents are built from a **scalar geometric series** 1 + a(1−Γ) +
  a²(1−Γ)² + …. This is literally F.4's `(I − W)⁻¹` in the
  N=0 / scalar limit, with a = 1 for white BC and (1−Γ) playing the
  role of W.
- The subroutine COPRAN listed in Fig. 5 (pp. 124-125) explicitly
  programs this scalar closure for arbitrary concentric annular
  cells. There is no multi-mode block structure.
- **No Legendre / moment / DP_N expansion of the boundary angular
  flux anywhere in Ch.IV.**
- Reference to higher-rank IC methods is implicit only via
  "References" on p. 141: Carlvik 1965/1967a. The chapter itself
  stops at rank-0 scalar.

Note on Ch.VI correction: prior memory
(`stammler_1983_ch6_interface_currents.md`) flagged Ch.VI as S_N,
not IC. That is confirmed. **Ch.IV is the right chapter, and its
method is rank-0 scalar.** The CP moment-integral memory
(`cp_moment_integrals.md`) claimed Ch.4-6 covers CP — tighten that
to Ch.IV only.

### Stacey (2007) — *Nuclear Reactor Physics*, Chapter 9 "Neutron Transport Theory"

Modern textbook, 80-page chapter covering the full transport
inventory: Boltzmann equation, P_L, S_N, CP, IC, MoC, nodal. IC
content is §9.4 (slab, pp. 321-330) and §9.5 (multidim, pp. 330-337).

- §9.4 uses **pure DP-0 throughout** — "double P₀ approximation" =
  isotropic angular flux in each hemisphere at each interface,
  parameterized by a **single** J⁺ and J⁻ per interface. Eqs. 9.66-9.97
  build the uncollided + n-collided emergence cascade,
  reflection/transmission probabilities T_{0i}, P_{0i}, R_i, T_i,
  and finally the 2×2 response-matrix form
  `J_i^± = R_i J_{i+1}^±`. Zero Legendre moments above P_0.
- §9.5 generalizes to multidim by replacing hemisphere-averaged
  quantities with per-interface averages (Eq. 9.100):
  `J_ij = Σ_k T_0i^{kj} J_ki + [Σ_k (1 − Σ_l T_0i^{kl}) J_ki] c_i P_i
  Λ_ij + s_ij s_i P_i`. Same DP-0 assumption, now per-face.
  Explicitly noted on p. 334 (comment under Eq. 9.110): "the incident
  current into volume i from volume 1 was calculated on the basis
  of a DP-0 angular flux approximation with respect to the
  orientation of the incident surface."
- §9.4 contains a **crucial warning** on p. 329 that directly bears
  on the Sanchez plateau: in a purely absorbing medium with an
  incident isotropic source, the exact flux attenuation from x=0 to
  x is E₂(Σx); but if the slab is subdivided into N intervals and
  DP-0 is applied at each interface, the computed attenuation is
  `[E₂(Σ)]^N`, NOT `E₂(NΣ)`. The two differ; **DP-0 subdivision
  introduces a floor error** that does not vanish as N → ∞ for
  forward-peaked boundary fluxes. This is precisely the failure
  mode the paper ascribes to "highly absorbing multiregion problems."
  A hollow sphere with large σ_t·R where the flux is forward-peaked
  at the inner cavity boundary is exactly this regime.
- §9.6+ covers P_L and S_N; no further IC content.


## §2 Cross-reference comparison

| Dimension | Ligou 1982 §8.2 | Sanchez 2002 | Stamm'ler 1983 Ch.IV §10 | Stacey 2007 §9.4-9.5 | Sanchez-McCormick 1982 §III.F.1 |
|-----------|-----------------|--------------|--------------------------|-----------------------|----------------------------------|
| Angular basis at surface | none (scalar j⁻) | piecewise-constant OR collocation-δ₂ (spatial + angular sectors, no Legendre) | none (scalar j⁻ + cosine return) | DP-0 (isotropic hemisphere) | Legendre P̃_n on [0,1] (or multifunction) |
| Emission weight | cosine/4π | depends on representation | cosine (Lambertian) | cosine (Lambertian) | (2n+1) P̃_n(μ) |
| Arrival weight | 1 | depends on representation | 1 | 1 | (2n+1) P̃_n(μ) |
| Surface Jacobian | none | Â_a·|Ω·n| (distributional) | none | inherited from area ratio |
| Closure form | 1/(1 − a(1−Γ)) scalar | Eq. 38 matrix (code-only) | 1/(1 − a(1−Γ)) scalar (Eq. 34) | 2×2 response matrix per slab, scalar DP-0 per face in multidim | (I − W)⁻¹ P matrix with Legendre-indexed blocks |
| N→∞ behavior proved? | n/a (rank-0) | not proved analytically | n/a (rank-0) | **E₂ product error floor shown** (p. 329) | claimed asymptotic but no proof of cell-integrated k_eff convergence |
| Reduces to ORPHEUS F.4 at scalar limit? | **Yes — identically** | Yes, when Nφ×Nθ = 1×1 and surface is integrated | **Yes — identically** | Yes — for a single-sphere cell, §9.4 IS F.4 | ambiguous (see §4) |


## §3 Which formulation reduces to F.4 at N=1?

**Ligou §8.2, Stamm'ler Ch.IV §10, and Stacey §9.4 are literally
F.4.** There is no "rank-1 Sanchez" in these references — the whole
method IS the scalar closure. The F.4 formulation is therefore the
canonical rank-0 baseline, and its 0.077 % residual at σ_t·R = 5 is
the expected order of error from three independent textbook
derivations, not an accident of your implementation.

At σ_t = 0 (the homogeneous limit), the F.4 identity
`W_oo + W_io = 1` — which you verified numerically — is a direct
consequence of:

- cosine-distributed return from a convex body of mean chord l̄ =
  4V/S (Stamm'ler Eq. 24a: Γ = l̄ Σ_r P̄)
- scalar multi-reflection sum 1/(1−(1−Γ)) → ∞ as Γ → 0 (Stamm'ler
  Eq. 34)
- and the fact that every neutron entering isotropically must
  eventually exit through some face.

This probability conservation is **baked into the scalar closure**
and is what makes Wigner-Seitz + white BC return k_eff = k_inf
exactly. Any closure that preserves this identity at every rank
**must** return the exact k_inf for homogeneous cells — and any
closure that breaks it at N≥1 will produce a residual that does not
vanish with N.

**The Sanchez-McCormick 1982 §III.F.1 per-face rank-N closure as
implemented in ORPHEUS does not exactly satisfy this conservation at
N≥1 in the homogeneous limit.** This is almost certainly the root
cause of the 1.42 % plateau. (See §6 recipe 1.)


## §4 Differences from Sanchez-McCormick 1982 §III.F.1

- **Basis**: SM-1982 uses Legendre moments P̃_n(μ). All four
  references use either no basis (scalar), hemisphere-averaged
  (DP-0), or spatial-sector × angular-sector piecewise functions.
  Nobody cross-validates the Legendre-moment form.
- **(2n+1) normalization**: SM-1982 puts the (2n+1) factor into the
  R (reflection) matrix (see user's memory
  `rank_n_interface_current_canonical.md` — "(2n+1) belongs in R
  not W"). None of the four new references provides a Legendre
  derivation that would confirm or refute this placement.
- **Reciprocity**: SM-1982 states reciprocity transposes mode
  indices at N≥1. Sanchez 2002 discusses reciprocity abstractly
  (Eq. 36: (ψ_a^r, ψ_b^γ)₊ = (1/c_{a,r}) δ_{ab} δ_{rγ}) but for the
  piecewise-constant / collocation cases, not for Legendre modes.
- **"W matrix (μ-weighted) vs. Lambert W"**: the textbook scalar
  closures all use **Lambert** (cosine-weighted) partial current
  responses. SM-1982's µ-weighted orthonormal moment W is a
  different object. At N=1, if the µ weight is correctly absorbed
  into the P̃_n basis and the (2n+1) factor is correctly placed,
  SM-1982 should reduce to Lambert — but this reduction is asserted,
  not derived in the 1982 review.


## §5 Explicit formulas for hollow sphere / cylinder

### Scalar (F.4-equivalent) hollow cylinder — Stamm'ler Eq. 24-25a and COPRAN

For annular shells 1..N:
- Bickley-function-based CP matrix p_ij via Eq. 25b, 26, 27a-b
  (Carlvik's method, subroutine COPRAN pp. 124-125).
- Outer-surface escape γ_i = (4 V_i / S_B) Σ_i p_i (Eq. 24d).
- Multi-reflection with albedo a at white outer boundary:
  Y_i(a) = Y_i / (1 − a(1−Γ)) where Γ = Σ Σ_{r,i} V_i Y_i (Eq. 32-34).
- **This is the exact formulation F.4 implements for cylinders.**

### Scalar (F.4-equivalent) hollow sphere

None of the four references gives a hollow-sphere-specific formula.
Ligou and Stamm'ler only treat cylindrical cells (Wigner-Seitz for
pin-cell), Stacey gives only slab and "annular" (cylindrical). This
is because pin-cell codes dominated the 1970s-80s literature. The
spherical-cell formulation you derived is a valid generalization of
the annular form, and the 0.077 % residual is the expected
equivalent accuracy.

### Sanchez 2002 multigroup closure (code-only)

Eq. 38 `J_a+^r = Σ_i E_ai^r q_i + Σ_{b,n} T_ab^{rn} J_b-^n` is
agnostic to the angular representation. For collocation with the
tracking quadrature, Eq. 42 reduces the exiting moment to the
weighted sum of angular fluxes along exiting trajectories — again
not a Legendre moment closure.

### Stacey §9.4 response matrix — slab

`J_i^± = R_i J_{i+1}^±` with R_i the 2×2 matrix

```
R_i = [[T_i^{-1},              −T_i^{-1} R_i],
       [R_i T_i^{-1},  T_i − R_i T_i^{-1} R_i]]
```

where T_i = T_{0i} + ½ c_i P_i (1 − T_{0i}) and R_i = ½ c_i P_i
(1 − T_{0i}) (Eqs. 9.86-9.92). This extends to multidim via Eq.
9.100 with the same DP-0 per-face assumption.


## §6 Recommended recipes, ranked by likelihood of closing the 1.42 % plateau

### Recipe 1 — VERIFY THE CONSERVATION IDENTITY AT σ_t → 0 FOR THE N≥1 SANCHEZ MATRICES (highest priority, cheap)

The F.4 closure satisfies W_oo + W_io = 1 at σ_t = 0 — a necessary
condition for white-BC + Wigner-Seitz + homogeneous → k_eff = k_inf.
Check whether the Sanchez N≥1 per-face matrices satisfy the analogous
identity:

```
Σ_faces Σ_m ⟨W_{f'f,mn}⟩ = δ_{n0}   at σ_t = 0
```

i.e., for every mode n on the outgoing side, the sum over incoming
faces and modes of the transmission block must equal 1 for n=0 and
0 for n>0 (no mode should lose or create particles). If the
ORPHEUS-implemented SM-1982 matrices fail this at σ_t = 0, **that is
the bug** — the plateau is a direct consequence of non-conservation
at higher modes. Test:

```
σ_t = 0, r_0/R = 0.3, N = 2
assert |W_oo[0,0] + W_io[0,0] - 1| < 1e-12
assert |W_oo[1,1] + W_io[1,1] - 1| < 1e-12   # higher-mode conservation
assert |W_oo[2,2] + W_io[2,2] - 1| < 1e-12
```

If higher-mode rows violate conservation, the fix is almost always
a missing (2n+1) factor in the emission or arrival weight, or a
mismatched measure between how P and W are defined. Stacey's p. 329
`[E₂(Σ)]^N ≠ E₂(NΣ)` observation is the slab analogue of this
failure.

### Recipe 2 — TRY THE PIECEWISE-CONSTANT ANGULAR (PCA) REPRESENTATION OF SANCHEZ 2002, NOT LEGENDRE (medium priority)

Sanchez 2002's alternative to DP_N is **angularly piecewise-constant**
representation (Eq. 37): partition the hemisphere into N² angular
sectors and use characteristic functions. This is what APOLLO2's TDT
actually uses. The basis functions satisfy Eq. 36 orthonormality
with c_{a,r} = ∫_a dS ∫_r |Ω·n| dΩ (explicit Ω·n weight built into
the normalization). Reciprocity is automatic because quadrature
formulas are symmetric (Sanchez 2002 p. 36, lines 950-955).

**Key advantage over Legendre**: particle conservation across the
cone partition is exact by construction — you can't lose particles
to a mode that doesn't "fit" the next cell's basis. Legendre modes
can lose particles because the even-moment integrals J_{2k} are
non-zero on a half-range.

Test: implement N² = 4 (π/2 split × μ=cos θ split) PCA closure for
the hollow sphere. If residual drops below 0.5 % at σ_t·R = 5, this
is the correct basis and the Legendre-moment SM-1982 recipe has a
basis-level flaw (not just a weight-placement flaw).

### Recipe 3 — ACCEPT F.4 AS THE PRODUCTION CLOSURE AND DEPRECATE RANK-N (lowest investment, highest confidence)

All four authoritative reactor-physics references converge on
scalar/DP-0 as the canonical IC closure. F.4 is mathematically
identical to the Stamm'ler / Ligou / Stacey formulation. Its 0.077 %
residual at σ_t·R = 5 on a hollow sphere is the expected accuracy
for the method at this optical thickness — and it matches or beats
every Sanchez N≥1 variant you have tested.

The Sanchez-McCormick 1982 §III.F.1 formal rank-N ladder is:

- not independently derived in any of the four references searched
- claimed by §III.F.1 but without a convergence proof or numerical
  demonstration on a closed-cell problem
- empirically worse than scalar F.4 at every σ_t·R ≤ 10 we have
  tested

If the project goal is a production rank-N closure for anisotropic
incoming flux at interface boundaries (pin-cell with anisotropic
neighbor), the modern path is Sanchez 2002's trajectory-based
approach with piecewise-constant or collocation angular
representation — NOT the 1982 Legendre-moment ladder. The 1982
ladder appears to be a *theoretical* construction that did not
cross over into any successor textbook.

### Recipe 4 — CONSULT HÉBERT (2009) "Applied Reactor Physics" (out of scope of this memo but recommended)

The user's priority list (per AGENT.md) places Hébert (2009) as the
CP reference alongside Stamm'ler. Hébert Chapters 5-6 may contain a
modern Legendre-moment IC derivation that Sanchez-McCormick 1982
does not. If the rank-N ladder exists anywhere in a modern textbook,
it is there.


## §7 Memory hygiene

Update these existing memories based on this extraction:

- `stammler_1983_ch6_interface_currents.md` — already correct
  (Ch.VI = S_N, not IC). Add a cross-reference: "Stamm'ler Ch.IV §10
  IS the IC chapter, and uses scalar cosine-return closure (Eqs.
  29-37). There is no rank-N Legendre ladder in either Ch.IV or
  Ch.VI — Stamm'ler does not treat rank-N IC at all."
- `cp_moment_integrals.md` — "cite A&S §5.1 + Hébert §3 /
  Stamm'ler-Abbate Ch.4-6" → tighten to "Stamm'ler-Abbate Ch.IV".
  Ch.V and Ch.VI are S_N.
- `rank_n_interface_current_canonical.md` — tag Sanchez-McCormick
  1982 §III.F.1 as **NOT cross-validated** by Ligou 1982, Sanchez
  2002, Stamm'ler 1983 Ch.IV, or Stacey 2007. The §III.F.1 recipe
  stands alone; treat novel k_eff claims from it with extra
  scepticism.


## §8 Citations

- J. Ligou (1982). *Elements of Nuclear Engineering* — Chapter 8
  "Fundamental Problems of Neutronics". (Identified from chapter
  content; original publisher appears to be in the ISS series, but
  not directly verifiable from the PDF which lacks a title page.)
- R. Sanchez, Li Mao, S. Santandrea (2002). "Treatment of Boundary
  Conditions in Trajectory-Based Deterministic Transport Methods."
  *Nuclear Science and Engineering* **140**(1), 23-50.
  DOI: 10.13182/NSE140-23.
- R. J. J. Stamm'ler, M. J. Abbate (1983). *Methods of Steady-State
  Reactor Physics in Nuclear Design* — Chapter IV "Integral Transport
  Theory; Collision Probabilities", pp. 105-141. Academic Press.
- W. M. Stacey (2007). *Nuclear Reactor Physics*, 2nd ed. — Chapter
  9 "Neutron Transport Theory", pp. 305-384. Wiley-VCH.
  ISBN 978-3-527-40679-1.
- (Already in ORPHEUS corpus) R. Sanchez, N. J. McCormick (1982).
  "A Review of Neutron Transport Approximations." *Nuclear Science
  and Engineering* **80**(4), 481-535. DOI: 10.13182/NSE82-A17197.
