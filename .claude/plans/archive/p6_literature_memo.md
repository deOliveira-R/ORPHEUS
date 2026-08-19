# P6 (#281) literature memo — adjoint-weighted (bilinear / eigenvalue-consistent) few-group-constant collapse

**Task**: extract the classical adjoint-weighted group-constant collapse
prescriptions (vector channels, scattering matrix two-flux rule, χ/fission
factored treatment, perturbation-theory grounding, condensation-vs-spatial
distinction) for ORPHEUS P6 adjoint-weighted homogenization/condensation.

**Method**: Tier-0 local library first (`scratch/literature/` + OCR sidecars
at `scratch/literature_ocr/`), per `.claude/rules/delegation.md`. Zotero MCP
tools are NOT exposed this session (no `mcp__zotero__*` in the tool surface)
— Tier 1 skipped, noted per L-007 discipline. Tier 2 (web) used ONLY for
provenance verification of citations found locally, never as a substitute
extraction source.

**Status**: COMPLETE (2026-07-25; EXTENDED 2026-07-26 with Bell &
Glasstone, acquired by the user into the local library — Source E).
Four local sources extracted + verified. The per-channel bilinear
prescriptions are now fully grounded in B&G Ch. 6 (§6.4h); remaining
acquisition candidates (Stacey Ch. 13, Williams CRC, Henry 1975) are
CORROBORATION-tier only, no longer blocking.

---

## 0. Local-library inventory verdict

Grep sweep of all sidecars for `adjoint`, `bilinear`, `perturbation`,
`condensation`, `homogeniz*`, `group constant`, `collapse`:

**FOUND LOCALLY (extractable content):**

| File | Relevant content |
|---|---|
| `Hebert(2009)Chapter3.pdf` | Adjoint transport equation, Rayleigh ratio + first-order stationarity statement, multigroup adjoint flux, flux-weighted condensation §3.5 (the baseline the bilinear collapse corrects) |
| `Nuclear Computational Science - A Century in Review.pdf` | Ch. 4 "Reactor Core Methods" §4.3.2 homogenization: importance-weighted reaction-rate conservation ⟨Φ†,Q⟩, flux-volume average, correction factors, GET/discontinuity factors; a later chapter (§8.5.4) explicitly states the asymptotic-homogenization result that homogenized constants are weighted by BOTH forward and adjoint flux |

**NOT IN LOCAL FOLDER (extraction STOPPED on each per the brief;
acquisition list at the end):**

- ~~Bell & Glasstone, *Nuclear Reactor Theory* (1970)~~ — **ACQUIRED
  2026-07-26** (`Bell-Glasstone(1970)Nuclear_reactor_theory.pdf`, 637-pp
  sidecar); extracted as Source E. The §6.4h bilinear machinery closed
  taxonomy rows 1, 2, 3, 4, 7.
- Stacey, *Nuclear Reactor Physics* — Ch. 13 "Perturbation and Variational
  Methods" (only Ch. 9 is local; its sole adjoint mention is a
  variance-reduction pointer to Ch. 13).
- Williams, "Perturbation Theory for Nuclear Reactor Analysis" (CRC
  Handbook of Nuclear Reactors Calculations, Vol. III, 1986).
- Duderstadt & Hamilton, *Nuclear Reactor Analysis* (1976) — Ch. 7
  few-group constants.
- Henry, *Nuclear-Reactor Analysis* (1975) — spectrum-weighted constants.
- Hébert (2009) chapters beyond Ch. 3 (the equivalence/SPH chapter is not
  local).
- Stammler & Abbate Ch. VI is local but (prior extraction, memory
  `energy-condensation-collapse-formulas`) contains only the flux-weighted
  prescription — re-checked this session, no adjoint weighting.

---

## Sources — extraction records

### Source A — Hébert (2009), *Applied Reactor Physics*, Ch. 3 [LOCAL: `Hebert(2009)Chapter3.pdf`]

Citation: A. Hébert, *Applied Reactor Physics*, Presses Internationales
Polytechnique, Montréal (2009), Chapter 3 "The Transport Equation".
All equations below spot-verified against the rendered scan (PDF pp. 11,
19–21 = printed pp. 77, 85–87).

**A.1 — The adjoint transport equation and the importance function
(§3.3.1, printed p. 77, PDF p. 11).**

Adjoint construction rules (paraphrase of the numbered list, p. 77):
(1) transpose matrix operators; (2) flip the sign of odd-parity
differential operators; (3) interchange the arguments of integral-operator
kernels.

- Eq. (3.60): −Ω·∇φ*(r,E,Ω) + Σ(r,E) φ*(r,E,Ω) = Q*(r,E,Ω).
- Eq. (3.61): the adjoint source density —

  Q*(r,E,Ω) = ∫₀^∞ dE′ Σ_ℓ (2ℓ+1)/(4π) Σ_{s,ℓ}(r, E′←E) Σ_m R_ℓ^m(Ω) φ*_ℓ^m(r,E′)
  + (1/(4π K_eff)) Σ_j νΣ_{f,j}(r,E) ∫₀^∞ dE′ χ_j(E′) φ*(r,E′).

  Two structural transpositions vs. the forward source (3.111-analogue):
  the scattering kernel keeps the FORWARD kernel Σ_{s,ℓ}(E′←E) but the
  free variable moves to the source slot (adjoint "scatters" E→E′ with
  the transposed kernel), and in the fission term νΣ_f sits at the LOCAL
  energy E while χ is integrated against φ* — the rank-1 fission dyad
  χ⊗νΣ_f transposes to νΣ_f⊗χ.
- Function-vs-distribution (p. 77, load-bearing for the χ channel):
  "the adjoint flux is a *function* of E" — it cannot be a distribution
  because χ_j(E′)φ*(r,E′) in (3.61) would otherwise multiply two
  distributions. "In some textbooks, the adjoint flux is referred to as
  the *importance function*." (short quotes, p. 77)
- Eq. (3.62): normalization of φ* is arbitrary, ∫dE ∫_V d³r φ* = 1.

**A.2 — Rayleigh ratio and first-order stationarity (§3.3.1, printed
p. 77, PDF p. 11) — the perturbation-theory grounding (brief item 4).**

- Eq. (3.63) (verified against the scan; it is written POINT-WISE in r,
  with NO spatial integral and NO leakage term — an infinite-medium /
  per-point form):

  K_eff = [ Σ_j ∫dE′ χ_j(E′) φ*(r,E′) · ∫dE νΣ_{f,j}(r,E) φ(r,E) ]
        / [ ∫dE′ φ*(r,E′) ( Σ(r,E′) φ(r,E′) − ∫dE Σ_{s,0}(r,E′←E) φ(r,E) ) ].

  Numerator = adjoint-weighted fission worth, factored as
  (∫χ_j φ* dE′)·(∫νΣ_{f,j} φ dE) per fissioning isotope j — the rank-1
  separability of the fission kernel appears explicitly.
  Denominator = adjoint-weighted (collision − in-scattering) = net
  removal worth.
- Stationarity statement (verbatim, p. 77, ≤3 lines): "The interest of
  the Rayleigh ratio is that it is stationary with respect to a small
  variation δΣ of the cross-section terms. The first order variation
  δK_eff … can be written in terms of φ, φ*, Σ and δΣ, *without* using
  δφ or δφ*. This observation is at the origin of the *perturbation
  theory*."
- What Hébert does NOT give in Ch. 3: the explicit δk formula
  δk ∝ ⟨φ*, δA φ⟩/⟨φ*, Fφ⟩ is asserted (existence) but not written out,
  and no bilinear-weighted condensation is derived from it.

**A.3 — Multigroup adjoint flux is a lethargy-AVERAGE, not an integral
(§3.5.1, printed p. 85, PDF p. 19).**

The forward multigroup flux is the group INTEGRAL (φ is a distribution);
the adjoint is defined as the group AVERAGE (φ* is a function):

- Eq. (3.118): φ*_g(r,Ω) ≡ ⟨φ*(r,Ω)⟩_g = (1/(u_g−u_{g−1})) ∫_{u_{g−1}}^{u_g} du φ*(r,u,Ω)
- Eq. (3.119): same for the angle-integrated φ*_g(r).

This is the same function/distribution split as Hébert's (3.96) [rate:
plain ∫_g du, no 1/Δu] vs (3.97) [function: lethargy average] — see
memory `energy-condensation-collapse-formulas`. For ORPHEUS: the adjoint
carrier condenses by AVERAGING, the forward by SUMMING. Duality pairing
⟨φ*,φ⟩ = Σ_g φ*_g φ_g is then exactly preserved when φ*_g is the
φ-weighted average over g (see caveat A.5).

**A.4 — The multigroup adjoint system = the TRANSPOSE of the multigroup
forward system, with the SAME (forward-flux-weighted) group constants
(§3.5.1, printed pp. 85–86, PDF pp. 19–20).**

- Eq. (3.120): −Ω·∇φ*_g + Σ_g(r) φ*_g = Q*_g.
- Eq. (3.121): Q*_g(r,Ω) = Σ_{h} Σ_ℓ (2ℓ+1)/(4π) Σ_{s,ℓ,h←g}(r) Σ_m R_ℓ^m(Ω) φ*_{ℓ,h}^m(r)
  + (1/(4π K_eff)) Σ_j νΣ_{f,j,g}(r) Σ_h χ_{j,h} φ*_h(r).

  Index transposition vs. the forward (3.111): scattering sums over SINK
  h with kernel Σ_{s,ℓ,h←g} (forward matrix, transposed indices);
  fission uses νΣ_f at the OWN group g times Σ_h χ_{j,h} φ*_h.
- CRITICAL GAP Hébert glosses over: (3.120)–(3.121) is "adjoint AFTER
  discretization" (transpose of the multigroup operator built with the
  §3.5 flux-weighted constants), yet Hébert *identifies* its solution
  with ⟨φ*⟩_g, the group-average of the continuous adjoint. That
  identification is exactly what bilinear-weighted constants are needed
  for; with flux-weighted constants "discretize-then-adjoint" ≠
  "adjoint-then-discretize" at first order. Hébert states no condition
  for the identification. (This is the classical B&G §6.4 discussion —
  source NOT local, see missing list.)

**A.5 — The flux-weighted baseline (§3.5.1) and the kinetics caveat
(§3.5.2, printed pp. 86–87, PDF pp. 20–21).**

Baseline (prior extraction, re-confirmed): (3.103) vector channels
Σ_g = ⟨Σφ⟩_g/φ_g; (3.104) scattering Σ_{s,ℓ,g←h} = ⟨Σ_{s,ℓ}φ⟩_{g←h}/φ_h
(SOURCE-flux weighted, sink-summed); (3.105) νΣ_{f,g} = ⟨νΣ_f φ⟩_g/φ_g;
(3.112) χ_{j,g} = ∫_g du χ_j(u) — plain sum, NO weight of any kind.
ALL forward-flux weighted; φ* appears in none of them.

Kinetics condensation (§3.5.2) is also purely flux-weighted:
- (3.125): 1/V_{n,g} = ⟨(1/V_n)φ⟩_g/φ_g — flux-weighted inverse velocity.
- (3.126)/(3.127): β_{ℓ,j} = Σ_g νΣ^del_{f,ℓ,j,g}φ_g / Σ_g νΣ_{f,j,g}φ_g.
- (3.131)–(3.133): isotope-merged β_ℓ, χ^pr_g, χ^del_{ℓ,g} — all weighted
  by the fission PRODUCTION rate νΣ_f φ, never by φ*.
- (3.135): Λ = Σ_g (1/V_{n,g})φ_g / Σ_j Σ_g νΣ_{f,j,g}φ_g.

Caveat for P6: these are the lattice-code conventions. The classical
*effective* kinetics parameters (β_eff, Λ_eff) are the adjoint-weighted
versions (bilinear ⟨φ*χ^del νΣ_f φ⟩-type ratios); Hébert Ch. 3 does not
introduce them. So Hébert supplies the adjoint OBJECTS (φ*, transposed
system, Rayleigh ratio) and the flux-weighted BASELINE, but never the
bilinear collapse itself.

---

---

### Source B — R. Roy, "Reactor Core Methods", Ch. 4 of Azmy & Sartori (eds.), *Nuclear Computational Science: A Century in Review*, Springer (2010) [LOCAL: `Nuclear Computational Science - A Century in Review.pdf`]

DOI of the volume chapters: 10.1007/978-90-481-3411-3 (chapter footers
carry per-chapter DOIs). §4.3.2 "Homogenization Process", printed
pp. 190–192, PDF pp. 200–202. Equations (4.41)–(4.46) spot-verified
against the rendered pages.

**B.1 — What the industrial prescription conserves (per Koebke [16]):**
volume-related conservation (integral flux + integral reaction rates in
the homogenized zone) and surface-related conservation (integral net
currents + integral fluxes at each interface). (paraphrase, p. 190)

**B.2 — Where the adjoint enters — as a RESPONSE weight, not a collapse
weight (§4.3.2.1, printed pp. 190–191):**

- Eq. (4.41): R_ρ = ⟨Σ_ρ, Φ⟩ = ∫_V d³r ∫_{4π} d²Ω Σ_ρ(r,Ω) Φ(r,Ω).
- Verbatim (p. 191, ≤2 lines): "the reaction rates are simply given by
  R_d = ⟨Φ†, Q⟩ with the appropriate importance function Φ†" — i.e.
  the lattice ADJOINT is invoked for detector-response functionals
  (duality ⟨Σ_d,Φ⟩ = ⟨Φ†,Q⟩), NOT as the weight inside the collapse.
- Eq. (4.42): conservation statement R_ρ = ⟨Σ_ρ,Φ⟩ = ⟨Σ̂_ρ, Φ̂⟩ (hat =
  homogeneous equivalent zone).
- Eq. (4.43): Σ̂_ρ ∫_V d³r Φ̂(r) = ∫_V d³r Σ_ρ(r)Φ(r) (isotropy assumed).
- Eq. (4.44): Σ̂_ρ = μ Σ̄_ρ = μ [∫Σ_ρΦ / ∫Φ] — flux-volume average TIMES
  a correction factor μ applied "to preserve reaction rates and integral
  values" when the homogeneous-equivalent flux differs from the
  reference (this μ is the SPH/equivalence-factor mechanism).
- Eqs. (4.45)–(4.46): net-current preservation at interfaces; Roy warns
  (p. 191, paraphrase) that net-flow-only conservation is degenerate in
  infinite-lattice calculations ("zero equals zero") and does not pin
  surface currents.
- Eq. (4.47) (§4.3.2.2, p. 192): GET discontinuity factors
  f^±_α = Φ_α,side/Φ̂_α,side per side of each interface [17 = K.S. Smith
  MIT PhD 1980].

**B.3 — Taxonomy verdict from Roy:** the production homogenization chain
is FORWARD-flux weighting + equivalence corrections (μ factors,
discontinuity factors), with the adjoint reserved for response
functionals. Bilinear collapse is absent from the industrial prescription
as Roy presents it. This is the "what the flux-weighted collapse gets
wrong is patched by equivalence factors, not by φ* weighting" branch of
the taxonomy.

---

### Source C — J. Dorning, "Nuclear Reactor Kinetics: 1934–1999 and Beyond", Ch. 8 of the same volume [LOCAL, same PDF]

Chapter DOI (printed footer, verified): 10.1007/978-90-481-3411-3_8.
Two distinct extractions: (C.1) the classical bilinear point-kinetics
parameters — the archetype of adjoint-weighted collapse-to-a-scalar; and
(C.2) the §8.5.4 statement that asymptotic spatial homogenization yields
FORWARD×ADJOINT-weighted constants. Equations (8.39)–(8.47) and the
§8.5.4 passages spot-verified against rendered PDF pp. 398–399, 443–444.

**C.1 — The Ussachoff/Henry adjoint-weighted point-kinetics parameters
(§8.3.2, printed pp. 388–392, PDF pp. 395–399).**

Lineage (paraphrase, p. 388): first reported "almost simultaneously" by
L.N. Ussachoff [27] and A.F. Henry [28], detailed derivation Henry [29].
Setup: n(r,v,Ω̂,t) = Ψ(r,v,Ω̂,t)·T(t) (Eq. 8.27, exact, not separation
of variables); reference-reactor forward/adjoint H₀N₀ = 0, H₀†N₀† = 0
(8.28/8.29); inner product ( , ) = ∫_V d³r ∫dv ∫d²Ω̂ (8.37);
normalization condition d/dt (N₀†, Ψ) = 0 (8.38). Projecting the
transport + precursor equations onto N₀† gives (8.35)/(8.36), and
division by the adjoint-weighted fission worth

- Eq. (8.39): F_t ≡ (N₀†, M_t Ψ) — THE bilinear normalization
  denominator ⟨φ*, Fφ⟩ of perturbation theory,
- with M_t = (1/4π) χ(v) ∫dv′∫d²Ω̂′ ν′... νΣ_f(r,v′,t) (8.40) and the
  steady-state spectrum χ = (1−β)χ₀ + Σᵢ βᵢχᵢ (8.41),

yields the point-kinetics equations (8.42)/(8.43) with the CLASSICAL
BILINEAR PARAMETERS:

- Eq. (8.44): ρ(t) = (1/F_t) (N₀†, (H_t − H₀) Ψ)
  — reactivity = ⟨φ*, δH φ⟩ / ⟨φ*, Fφ⟩. This IS the first-order
  perturbation δk formula of brief item 4, obtained here exactly (with
  Ψ the instantaneous shape) by subtracting (N₀†, H₀Ψ) =
  (H₀†N₀†, Ψ) = 0 (p. 391 — the adjoint kills the unperturbed term).
- Eq. (8.45): β̄ᵢᵗ = (1/F_t)(N₀†, M_tⁱ Ψ) — effective delayed fraction:
  the delayed-emission operator M_tⁱ (spectrum χᵢ) bilinearly weighted
  and normalized by the TOTAL fission worth.
- Eq. (8.46): ℓ̄_t = (1/F_t)(N₀†, Ψ) — effective generation time (in
  number-density form; in flux variables this is the ⟨φ*, (1/v)φ⟩ /
  ⟨φ*, Fφ⟩ inverse-velocity bilinear).
- Eq. (8.47): S̄(t) = (N₀†, S)/(N₀†, Ψ).

Interpretation (paraphrase p. 391): N₀† is "the neutron importance";
T(t) ∝ (N₀†, n) = total importance carried by the instantaneous
population. NOTE for the χ channel: in M_t the fission dyad is
χ(v)⊗νΣ_f(v′); its bilinear contraction (N₀†, M_tΨ) factorizes per
point as [∫dv N₀†χ]·[∫dv′ νΣ_fΨ] — adjoint contracts the EMISSION
spectrum, forward contracts the PRODUCTION cross-section. Same
factored structure as Hébert (3.63) numerator (Source A.2).

Dorning's caveat (paraphrase, p. 392): the normalization condition
(8.38) restricts the shape function in a "not-very-well-defined way";
variational [39 = M. Becker, "A generalized formulation of point nuclear
reactor kinetics equations", NSE 31(3):458–463 (1968); cf. [37] =
J. Lewins, J. Nucl. Energy A 12:108 (1960)] and asymptotic [40, 42 =
Dorning-Spiga 1978; Dorning 1980, ANS Transactions] formulations are
"truer approximations". NOTE: the whole chapter contains NO explicit
"errors are second order / stationary" sentence (chapter-wide grep for
{second order, stationary, insensitive} — zero relevant hits) — the
first-order-stationarity grounding in the local corpus is therefore
Hébert's p. 77 statement (A.2), not Dorning.

**C.2 — §8.5.4 "Homogenization Theories …" (printed pp. 435–437, PDF
pp. 442–444): asymptotic homogenization PRODUCES bilinear weighting.**

Two "approaches to homogenization" for nodal methods (paraphrase):

1. Exact/generalized equivalence theory [74 = K.S. Smith, Prog. Nucl.
   Energy 17(3):303–335 (1986); 76 = Koebke, IAEA TCM Vienna (1978);
   103 = Koebke (1981)]: homogenized constants = node-interior
   FORWARD-flux-weighted averages + discontinuity factors carry the
   equivalence burden. (Matches Roy, Source B.)
2. Multiscale asymptotic expansion — THE key verbatim finding (printed
   p. 436, PDF p. 443, ≤3 lines):

   > "…the expressions for the assembly-homogenized cross sections and
   > diffusion coefficient are not given by the transport-flux-weighted
   > averages … rather, they are weighted by both the transport theory
   > flux and the transport theory adjoint flux (or importance) — a
   > result that follows directly from the solvability conditions"

   [refs 79, 80]. The bilinear weight is NOT a modeling choice here —
   it FALLS OUT of the Fredholm-alternative solvability conditions of
   the two-scale expansion.

3. The diffusion-based variant [106] replaces transport forward/adjoint
   by diffusion forward/adjoint in the same expressions (p. 436,
   paraphrase).
4. The three-scale theory [81] (printed p. 437, PDF p. 444, paraphrase):
   cell-level constants = cell averages "weighted by transport theory
   forward and adjoint flux solutions for the heterogeneous cell";
   assembly-level constants = assembly averages weighted by
   cell-homogenized diffusion FORWARD and ADJOINT solutions; each level
   also gets discontinuity factors. Fully recursive bilinear weighting.

**C.3 — Primary sources named by Dorning for the bilinear-weighting
result (chapter-8 reference list, sidecar-extracted, identities below;
NONE in the local folder — acquisition candidates):**

- [77] E.W. Larsen, "Neutron transport and diffusion in inhomogeneous
  media. I", J. Math. Phys. 16(7):1421–1427 (1975).
  DOI 10.1063/1.522714 [CrossRef-verified].
- [78] E.W. Larsen, "Neutron Transport and Diffusion in Inhomogeneous
  Media. II", Nucl. Sci. Eng. 60(4):357–368 (1976).
  DOI 10.13182/NSE76-A26897 [CrossRef-verified].
- [79] R.T. Chiang & J.J. Dorning, "A homogenization theory for lattices
  with burnup and non-uniform loadings", Proc. ANS 1980 Adv. React.
  Phys. Shielding, pp. 319–329 (1980). [conference; not DOI-resolvable]
- [80] H. Zhang, Rizwan-uddin, J.J. Dorning, "Systematic homogenization
  and self-consistent flux and pin power reconstruction for nodal
  diffusion methods — Part II", Transp. Theory Stat. Phys. 26(4):433–468
  (1997). DOI 10.1080/00411459708017925 [CrossRef-verified].
- [81] H. Zhang, Rizwan-uddin, J.J. Dorning, "A multiple-scales
  systematic theory for the simultaneous homogenization of lattice cells
  and fuel assemblies", Transp. Theory Stat. Phys. 26(7):765–811 (1997).
  DOI 10.1080/00411459708224422 [CrossRef-verified].
- [106] H. Zhang, Rizwan-uddin, J.J. Dorning, "Systematic Homogenization
  and Self-Consistent Flux and Pin Power Reconstruction for Nodal
  Diffusion Methods — 1: Diffusion Equation-Based Theory", Nucl. Sci.
  Eng. 121(2):226–244 (1995). DOI 10.13182/NSE95-A28560
  [CrossRef-verified].
- [74] K.S. Smith, "Assembly homogenization techniques for light water
  reactor analysis", Prog. Nucl. Energy 17(3):303–335 (1986).
  DOI 10.1016/0149-1970(86)90035-1 [CrossRef-verified].
- Kinetics lineage: [27] L.N. Ussachoff, "Equation for the importance of
  neutrons, kinetics and the theory of perturbations", Proc. Int. Conf.
  Peaceful Uses of Atomic Energy, Geneva 1955, P/656, Vol. 5,
  pp. 503–510; [28] A.F. Henry, WAPD-124 (1955); [29] A.F. Henry,
  "The Application of Reactor Kinetics to the Analysis of Experiments",
  Nucl. Sci. Eng. 3:52–70 (1958). DOI 10.13182/NSE58-1
  [CrossRef-verified].

---

### Source D — negative results (full-corpus sweep record)

Every sidecar was grep-swept for {adjoint, bilinear, perturbation,
condensation, homogeniz*, group constant, collapse, importance}. Files
with hits that are IRRELEVANT to the bilinear collapse (verified by
reading context): Adams-Larsen 2002 (iteration theory; "adjoint" =
transport-sweep machinery), Sanchez 1982 review (single hit = a cited
title), Stacey Ch. 9 (single hit = Monte-Carlo variance-reduction
pointer to its Ch. 13 — NOT local), Stammler Ch. IV/VI (zero adjoint
content; Ch. VI collapse is purely flux-weighted, confirming memory
`energy-condensation-collapse-formulas`), Brockmann 1981, Benoist 1981
(Wigner-Seitz diffusion coefficients — flux/current-weighted, not
adjoint), Hammer-Morel-Wang 2019, Valougeorgis 1985, Mitsis 1963,
Alcouffe 1977, Morel 1989, Sood 1999/2003, Larsen-Morel-Miller 1987,
Pomraning 1989. Roy's §4.3.1 (4.9)–(4.10) is Bondarenko
flux/(σ_t+σ_0)-weighted self-shielding — not adjoint. No NSE paper on
generalized equivalence / SPH / adjoint-weighted condensation is
currently in the folder.

---

### Source E — Bell & Glasstone (1970), *Nuclear Reactor Theory*, Ch. 6 [LOCAL since 2026-07-26: `Bell-Glasstone(1970)Nuclear_reactor_theory.pdf`] — the per-channel bilinear prescriptions

Citation: G.I. Bell & S. Glasstone, *Nuclear Reactor Theory*, Van
Nostrand Reinhold, New York (1970), TID-25606. Sidecar printed-page
mapping: printed ≈ PDF − 18. ALL load-bearing equations below
spot-verified against rendered PDF pp. 290, 297, 323–326 (printed
pp. 272, 279, 305–308). This section CLOSES the acquisition-blocked
items of the 2026-07-25 memo (taxonomy rows 1, 2, 3, 7).

**E.0 — Where it lives.** The bilinear collapse is §6.4h "Variational
Derivation of Multigroup Equations" (printed pp. 305–308) + §6.4i
"Self-Consistent Determination of Group Constants" (pp. 308–310),
grounded by §6.4b (the second-order functional theorem, pp. 292–295),
§6.4c (stationary eigenvalue functionals, p. 295), §6.2c (the
consistency failure of flux-weighted constants, pp. 272–273), and §6.3c
(the δk formula, pp. 277–279). Setting: P₁, plane geometry, fixed-source
functional J of Eq. (6.120); within-group separability is the ANSATZ.

**E.1 — The trial-function ansatz and the carrier normalizations
(§6.4h, printed p. 305, PDF p. 323).**

Trial functions (verified):
- (6.123): Φ(x,μ,E) ≃ φ_{0,g}(x)ψ_{0,g}(E) + 3μ φ_{1,g}(x)ψ_{1,g}(E)
- (6.124): Φ†(x,μ,E) ≃ φ†_{0,g}(x)ψ†_{0,g}(E) + 3μ φ†_{1,g}(x)ψ†_{1,g}(E)

with E_g ≤ E ≤ E_{g−1}; ψ, ψ† = KNOWN within-group spectra (e.g. from
infinite-medium or B₁ calculations, §4.5c). Normalizations (verified):

- (6.125): ∫_g ψ_{0,g}(E) dE = 1        ← the FORWARD shape integrates to 1
- (6.126): ∫_g ψ_{0,g}(E) ψ†_{0,g}(E) dE = 1   ← unit BILINEAR overlap
- (6.127): ∫_g ψ_{0,g}(E) ψ†_{1,g}(E) dE = 1
- (6.128): ∫_g ψ†_{0,g}(E) ψ_{1,g}(E) dE = 1

Verbatim (printed p. 306, PDF p. 324, ≤2 lines): "From equation
(6.125), it is seen that by this normalization φ_{0,g}(x) is simply
1/4π times the total flux of neutrons in group g." Flat-spectrum limit
(same page): ψ_{0,g} ≃ ψ_{1,g} ≃ 1/ΔE_g and ψ†_{0,g} ≃ ψ†_{1,g} ≃ 1.

**THE CARRIER VERDICT (brief question 2).** B&G's coarse flux carrier
is the PLAIN condensed flux: φ_{0,g} = (1/4π)·∫_g φ dE — a group
INTEGRAL, exactly Hébert's forward convention. The coarse adjoint
carrier φ†_{0,g} is scaled by the unit-bilinear-overlap condition
(6.126). Synthesis (arithmetic from their stated normalizations, not
their prose): since ⟨ψ†_{0,g}⟩_{ψ₀-weighted} = ∫ψ†₀ψ₀/∫ψ₀ = 1/1 = 1,
the coarse adjoint carrier is the **FLUX-WEIGHTED group-average of the
pointwise adjoint**, Φ*_G = ⟨φ*φ⟩_G/Φ_G. Consequences:

- The duality pairing is exactly preserved by construction:
  ⟨φ*φ⟩_G = Φ*_G · Φ_G (per group, per moment, via (6.126)–(6.128)).
- Under the within-group separability ansatz, ORPHEUS conventions (a)
  [plain-Φ carrier, Σ_C = ⟨φ*Σφ⟩_G/(Φ*_G·Φ_G)] and (b) [diagonal
  bilinear ⟨φ*Σφ⟩_G/⟨φ*φ⟩_G with bilinear carrier] COINCIDE
  numerically — B&G's normalization makes the (a)-denominator equal
  the (b)-denominator. B&G's ROW EQUATIONS are written in the (a)
  shape: plain carriers φ_{0,g}, φ_{1,g} multiplied by bilinear-valued
  constants (see (6.131)/(6.132) below). The (a)/(b) distinction only
  reopens when within-group shape is NOT separable in (x,E) — which is
  exactly what §6.4i's region-wise self-consistent spectra address.
- Contrast with Hébert (3.118) (Source A.3): Hébert's multigroup
  adjoint is the PLAIN lethargy-average; B&G's is the flux-weighted
  average. They agree in the flat-within-group-flux limit. B&G's is
  the dual-consistent one (preserves ⟨φ*,φ⟩ exactly).

**E.2 — The vector-channel bilinear rule (brief question 1). Eq.
(6.135), printed p. 307, PDF p. 325 (verified):**

  σ_{i,g}(x) ≡ ∫_g σ(x,E) ψ†_{i,g}(E) ψ_{i,g}(E) dE   (i = 0, 1)

- YES, this is the diagonal-pair bilinear: with normalization (6.126)
  the denominator ⟨ψ†ψ⟩_g = 1 is ABSORBED, so numerically
  σ_{0,g} = ⟨ψ†σψ⟩_g/⟨ψ†ψ⟩_g — the ⟨φ*Σφ⟩_G/⟨φ*φ⟩_G form.
- MOMENT-RESOLVED: i = 0 (flux channel, weight ψ†₀ψ₀) and i = 1
  (current channel, weight ψ†₁ψ₁) give DIFFERENT collapsed totals —
  the classical moment-dependent σ_g pair (the P₁ transport-correction
  axis), one per Legendre row.

**E.3 — The scattering-matrix two-flux rule (brief question 3). Eq.
(6.136), printed p. 307, PDF p. 325 (verified):**

  σ_{i,g′→g}(x) ≡ ∫_g dE ∫_{g′} dE′ σ_i(x; E′→E) ψ_{i,g′}(E′) ψ†_{i,g}(E)   (i = 0, 1)

YES — per-(g′,g)-pair: SOURCE-group forward spectrum ψ_{i,g′}(E′) ×
SINK-group adjoint spectrum ψ†_{i,g}(E). No denominator (absorbed by
(6.126)–(6.128)). This adjudicates the P6 load-bearing question:
the classical prescription IS sink-adjoint × source-flux per pair,
for every Legendre moment. (Baseline contrast: Ch. 4's (4.27) is
source-flux-only, ψ† ≡ 1.)

The resulting multigroup system (verified, printed p. 307):
- (6.131): d/dx[φ_{1,g}] + σ_{0,g}φ_{0,g} = Σ_{g′} σ_{0,g′→g}φ_{0,g′} + Q_{0,g}
- (6.132): d/dx[φ_{0,g}] + 3σ_{1,g}φ_{1,g} = 3Σ_{g′} σ_{1,g′→g}φ_{1,g′} + 3Q_{1,g}
- Sources are ADJOINT-spectrum-weighted: (6.133)/(6.134)
  Q_{i,g} = ½∫_g ψ†_{i,g}(E) dE ∫ P_i(μ) Q dμ — the external source
  collapses with the IMPORTANCE weight, not the flux weight.
- Verbatim (p. 307, ≤2 lines): "equations (6.131) and (6.132) are
  identical in form with the multigroup P₁ equations, but the group
  constants are now defined with both flux and adjoint (importance)
  weighting."
- The adjoint system (6.137)/(6.138) (printed pp. 307–308) uses THE
  SAME constants with g↔g′ transposed, and "equations (6.137) and
  (6.138) are clearly adjoint to equations (6.131) and (6.132)"
  (p. 308) — with bilinear constants, discretize-then-adjoint =
  adjoint-then-discretize. Structural reading (my language): §6.4h is
  a Petrov-Galerkin discretization in energy — trial basis = ψ
  (within-group flux spectra), test basis = ψ† (within-group adjoint
  spectra); (6.131)/(6.137) are the test-/trial-side projections of
  one bilinear form.

**E.4 — The fission/χ treatment (brief question 4).**

B&G carry fission INSIDE the transfer kernel: Ch. 4, Eq. (4.38)
(printed p. 186, sidecar-verified): σ_{0,g′→g} = σ_{s0,g′→g} +
νσ_{f,g′→g}, "the energy spectrum of these neutrons can be described
as part of the term σ_{0,g′→g}" (p. 186); the k-eigenvalue system
(4.39)/(4.41) divides only νσ_{f,g′→g} by k. Therefore under §6.4h the
fission channel collapses per-pair EXACTLY like scattering, by (6.136)
with σ_0 → νσ_f(E′→E):

  νσ_{f,g′→g} = ∫_g dE ∫_{g′} dE′ νσ_f(x; E′→E) ψ_{0,g′}(E′) ψ†_{0,g}(E)

Synthesis (grounded in (6.136)+(4.38), not spelled out by B&G): for the
separable kernel νσ_f(E′→E) = χ(E)·νσ_f(E′) this FACTORIZES per pair,

  νσ_{f,g′→g} = [∫_g χ(E)ψ†_{0,g}(E) dE] · [∫_{g′} νσ_f(E′)ψ_{0,g′}(E′) dE′]
              = χ†_g · (νσ_f)_{g′},

i.e. the collapsed spectrum is the ADJOINT-CONTRACTED (importance-
weighted) emission fraction χ†_g and the collapsed production is the
FORWARD-flux-weighted (νσ_f)_{g′} — the rank-1 dyad STAYS FACTORED,
sink-adjoint ⊗ source-flux. Consistently, the δk normalization (6.71)
(E.5) is the bilinear fission worth ⟨Φ†, (1/4π)νσ_f(E′→E) Φ⟩. Note
Σ_g χ†_g ≠ 1 in general — the adjoint-contracted spectrum is NOT a
unit-normalized distribution; its magnitude carries fission importance
(cf. Hébert (3.112) plain-sum χ, which keeps Σχ = 1 but is first-order
k-inconsistent).

**E.5 — Perturbation grounding and consistency statements (brief
questions 4-of-original + 5).**

- δk formula, §6.3c Eq. (6.71) (printed p. 279, PDF p. 297, verified):

  (Δk/k*) ∬∬ (1/4π)νσ_f(r,E′→E) Φ†(r,Ω,E) Φ(r,Ω′,E′) dV dΩ′dE′ dΩ dE
    ≈ −∭ Δσ(r,E) Φ†ΦdV dΩ dE + ∬∬ Δ[σf(r;Ω′,E′→Ω,E)] Φ†(Ω,E)Φ(Ω′,E′) …

  = the exact ⟨φ*, δA φ⟩/⟨φ*, Fφ⟩ structure (first order: Φ* → Φ;
  k* ≈ 1). Note the transfer-perturbation term carries the SAME
  sink-adjoint × source-flux pairing as (6.136).
- The second-order theorem, §6.4b Eq. (6.90) (printed p. 293, PDF
  p. 311, sidecar-verified): J = J₀ + (δΦ†, L δΦ). Verbatim (≤3
  lines): "the correction would be second order in small quantities …
  a better estimate than could be obtained, for example, from the
  inner product (Q†, Φ), since this would have an error (Q†, δΦ),
  which is first order in small quantities." THIS is the formal
  statement that bilinear-weighted (stationary-functional) collapse
  zeroes the first-order spectrum-error, while flux-only weighting
  retains it.
- Eigenvalue version, §6.4c Eq. (6.92): α ≃ (Φ†, LΦ)/[(1/v)Φ†, Φ] —
  the stationary Rayleigh-quotient estimate (k analog with the fission
  operator; used with trial spectra it makes the collapsed eigenvalue
  second-order accurate in δψ, δψ†).
- Dual consistency: §6.4g closing sentence (printed p. 305, verified):
  the variational route guarantees "the equation containing the
  approximate Φ† will indeed be adjoint to that for the approximate
  Φ, a result which is frequently not true for the equations derived
  in a more straightforward way." And §6.2c (printed p. 272, PDF
  p. 290, verified, the statement that closes Source A.4's gap):
  "the adjoint equations (6.56) and (6.57) could not have been derived
  by integrating an energy-dependent adjoint equation over an energy
  interval … the required cross sections, which are flux-weighted
  averages, would not be obtained. This problem is examined in §6.4h."
  (Same page, historical: "adjoint difference equations" of certain
  S_N codes are "not quite adjoint" to the flux difference equations
  in curved geometry — the discrete-adjoint hazard is old.)
- Practical status + trade-off (printed p. 308, PDF p. 326, verified,
  the balance caveat): "Relatively little use has yet been made of
  bilinearly (flux and adjoint) averaged group constants, primarily
  because the group adjoints must be estimated in addition to the
  group fluxes. When the bilinear averaging has been used, however,
  it seems to be superior to simple flux averaging, at least for
  problems involving only a few groups.³⁵ When a large number of
  groups can be employed, the adjoint weighting is less important,
  because the adjoint function will not vary significantly across a
  group."
- Rate-preservation trade-off (structural synthesis, flagged as such):
  (6.131) is the ψ†-weighted projection of the balance equation, so
  with bilinear σ_{0,g} the product σ_{0,g}φ_{0,g} is an
  importance-weighted removal, NOT the true collision rate (the true
  rate needs the Ch.-4 flux-weighted σ). B&G keep the carrier
  rate-meaningful ((6.125): φ_{0,g} IS the true group flux) and put
  ALL the bilinear content in the constants; the ψ† → 1 limit
  recovers Ch. 4's flux-weighted constants and the true balance. The
  k-functional is what the bilinear system preserves to second order —
  not channel-wise reaction rates.
- B&G's OWN kinetics chapter consumes §6.4h: Ch. 9 derives point
  kinetics with "bilinear weighted cross sections (§6.4h)" (sidecar
  line 14648) — the Henry lineage is internally consistent.
- Caveats B&G state about the variational estimate (printed
  pp. 293–294): the sign of (δΦ†, LδΦ) is unknown for
  energy-dependent problems (no minimum principle, unlike one-speed
  isotropic); L is unbounded so small δΦ does not imply small LδΦ;
  J depends on trial normalization (fix: the Schwinger form J_s,
  Eq. (6.91)).

**E.6 — Ch. 4 flux-weighted baseline (for the taxonomy contrast;
sidecar-verified):** (4.26) σ_{n,g} = ∫_g σφ_n dE / φ_{n,g};
(4.27) σ_{n,g′→g} = ∫_g dE ∫_{g′} dE′ σ_n(E′→E)φ_n(E′) / φ_{n,g′} —
source-moment-weighted only; (4.38) fission-in-transfer split. B&G
already weight per Legendre MOMENT n even in the flux-weighted
baseline (φ_0 = flux, φ_1 = current — Stammler's VI(6d) current-
weighted P1 channel matches this).

**E.7 — Reference identities (Ch. 6 list, sidecar lines 10055–10065):**
- [34] A.F. Henry, Nucl. Sci. Eng. 27:493 (1967) — time-dependent
  (delayed-neutron) generalization of the bilinear multigroup.
- [35] T.A. Pitterle & C.W. Maynard, Trans. ANS 8(1):205 (1965);
  W.W. Little Jr. & R.W. Hardie, NSE 29:402 (1967); A.J. Buslik, NSE
  32:233 (1968); "see, however" J.B. Yasinsky & S. Kaplan, NSE 31:354
  (1968) — the empirical basis (+ the flagged counterpoint) for the
  few-group-superiority claim.
- [36] Pomraning & Clark (their Ref. 30) — boundary trial-function
  restriction (6.139a).
- [37] T. Toivanen, J. Nucl. Energy 22:283 (1968); Lancefield (Ref.
  31); A.J. Buslik, Trans. ANS 12:152 (1969); R.B. Nicholson, Trans.
  ANS 12:731 (1969); R.J. Neuhold & K.O. Ott, NSE 39:14 (1970) —
  self-consistent within-group-spectra applications (§6.4i).

---

## Consolidated taxonomy table

Channel class → classical prescription → where it is grounded LOCALLY →
status. *(Updated 2026-07-26: B&G acquired and extracted — Source E.
Rows 1, 2, 3, 4, 7 last-column cells are now fully LOCAL.)*

| # | Channel class | Classical (bilinear / eigenvalue-consistent) prescription | Local grounding (verified) | Per-channel formula w/ eq. number |
|---|---|---|---|---|
| 1 | Vector channels (Σ_t, Σ_a, Σ_f, νΣ_f) | Preserve the adjoint-weighted reaction worth: Σ_G from ⟨φ*Σφ⟩_G/⟨φ*φ⟩_G-type ratios so that δk of the collapse vanishes at first order | Stationarity of K_eff w.r.t. δΣ at fixed φ,φ*: Hébert (3.63)+p. 77 statement. "Weighted by both the transport theory flux and the … adjoint flux" for homogenized constants: Dorning §8.5.4 p. 436 (from Fredholm solvability, refs [77]–[81]) | **LOCAL — B&G (6.135)**, p. 307: σ_{i,g} = ∫_g σ ψ†_{i,g} ψ_{i,g} dE (i=0,1; per Legendre moment) ≡ ⟨φ*Σφ⟩_G/⟨φ*φ⟩_G via the unit-overlap normalization (6.126). Flux-weighted baseline (3.103)/B&G (4.26) |
| 2 | Scattering matrix Σ_s(g′→g) | Two-flux rule: each (g′,g) pair weighted by SOURCE forward flux φ_{g′} AND SINK adjoint φ*_g (the pairing in ⟨φ*,Sφ⟩ = Σ_g Σ_{g′} φ*_g Σ_{s,g←g′} φ_{g′}) | The bilinear pairing structure is explicit in Hébert (3.61) continuous / (3.121) multigroup (transposed kernel: sink-adjoint contraction) and Dorning (8.30)+(8.37). No local source writes the fine→coarse two-flux collapse | **LOCAL — B&G (6.136)**, p. 307: σ_{i,g′→g} = ∫_g dE ∫_{g′} dE′ σ_i(E′→E) ψ_{i,g′}(E′) ψ†_{i,g}(E) — CONFIRMED per-pair SINK-adjoint × SOURCE-flux, both moments. Baseline (3.104)/B&G (4.27): source-flux only |
| 3 | Fission source: χ and νΣ_f | Factored rank-1 treatment: adjoint contracts the EMISSION spectrum (⟨χφ*⟩), forward contracts PRODUCTION (⟨νΣ_f φ⟩); preserved invariant = total fission worth ⟨φ*,Fφ⟩ | Rank-1 factorization is explicit in Hébert (3.61), (3.63) numerator = Σ_j ⟨χ_j φ*⟩·⟨νΣ_{f,j}φ⟩, and Dorning F_t (8.39)–(8.41); worked collapsed-spectrum-channel example = β̄ᵢ (8.45) = (N₀†,M_tⁱΨ)/F_t | **LOCAL — B&G (4.38)+(6.136)**: fission is a TRANSFER, collapses per-pair; separable kernel ⟹ νσ_{f,g′→g} = χ†_g·(νσ_f)_{g′}, χ†_g = ∫_g χ ψ†_{0,g} dE (adjoint-contracted; Σ_g χ†_g ≠ 1), (νσ_f)_{g′} flux-weighted — dyad stays factored (E.4). Flux-weighted baseline: χ_{j,g} plain sum (3.112) |
| 4 | Perturbation-theory grounding (δk formula) | δρ = ⟨φ*, δH φ⟩/⟨φ*, Fφ⟩; bilinear-weighted constants zero the first-order collapse error because K_eff is stationary in (φ,φ*) | Dorning (8.44) ρ = (N₀†,(H_t−H₀)Ψ)/F_t with F_t = (N₀†,M_tΨ) — the formula EXACTLY, reactivity flavor; Hébert (3.63) + the p. 77 stationarity statement ("δK_eff … without using δφ or δφ*") | **LOCAL (all three)**. δk formula: B&G (6.71) p. 279. The second-order THEOREM: B&G (6.90) p. 293, J = J₀ + (δΦ†, LδΦ) — flux-only weighting errs at FIRST order, (Q†, δΦ). Stationary eigenvalue functional: (6.92) |
| 5 | Condensation (energy) vs homogenization (space) | Energy condensation: classical texts give flux-weighted collapse + bilinear correction from perturbation theory. Spatial homogenization: EITHER flux-weight + equivalence factors (μ/SPH, discontinuity factors) OR asymptotic multiscale which FORCES forward×adjoint weighting | Roy §4.3.2 (4.41)–(4.47): industrial spatial chain = flux-volume average (4.44) + μ factors + GET f^± (4.47); adjoint only for response functionals ⟨Φ†,Q⟩. Dorning §8.5.4: asymptotic route yields bilinear weights + discontinuity factors at every scale (2- and 3-scale) | Hébert §3.5 = the flux-weighted energy-condensation baseline (LOCAL). What flux weighting gets wrong: **B&G p. 308** — bilinear "superior … only a few groups"; many groups ⟹ φ* ≈ flat within a group and the correction fades [refs 35, incl. Yasinsky-Kaplan counterpoint]; first-order error mechanism = (6.90) |
| 6 | (bonus) Kinetics parameters β_eff, Λ | β̄ᵢ = ⟨φ*,M^iφ⟩/⟨φ*,Fφ⟩, Λ̄ = ⟨φ*,(1/v)φ⟩/⟨φ*,Fφ⟩ — the oldest adjoint-weighted collapse (to a 0-D model) | Dorning (8.45)/(8.46) [LOCAL, verified]; lineage Ussachoff 1955 / Henry 1955, 1958 | LOCAL in operator form. Hébert's lattice-code (3.126)–(3.135) are the flux-weighted counterparts — a live example of the two conventions coexisting |
| 7 | (bonus) Multigroup-adjoint consistency | "Transpose the multigroup operator" ≠ "discretize the continuous adjoint" unless constants are consistently weighted; adjoint carrier is a group-AVERAGE (function), forward a group-INTEGRAL (distribution) | Hébert (3.118)–(3.121) + the distribution/function statements pp. 77, 85 [LOCAL, verified]. Hébert does NOT state the consistency condition | **LOCAL — B&G §6.2c p. 272** (verbatim in E.5): flux-weighted multigroup adjoint "could not have been derived by integrating an energy-dependent adjoint equation"; §6.4h fixes it — (6.137)/(6.138) share the (6.135)/(6.136) constants and are "clearly adjoint" (p. 308); §6.4g: the variational route guarantees dual consistency |

**Agreement/disagreement between the local sources**: no conflicts —
the four local sources partition cleanly: Hébert = adjoint objects +
flux-weighted collapse baseline; Roy = industrial spatial homogenization
(flux weight + equivalence corrections, adjoint for responses only);
Dorning = bilinear kinetics collapse + spatial homogenization as a
solvability RESULT; **B&G = the per-channel bilinear collapse
prescriptions themselves** ((6.135)/(6.136) + carriers (6.125)–(6.128)
+ theorem (6.90)). The load-bearing scattering question IS now
adjudicated: **B&G (6.136) prescribes per-pair sink-adjoint ×
source-flux weighting**, for every Legendre moment. Carrier
adjudication: B&G write the row equations in PLAIN carriers (forward =
group integral, (6.125)) with the coarse adjoint = FLUX-WEIGHTED
group-average (from (6.126)), which makes ⟨φ*φ⟩_G = Φ*_G·Φ_G exactly —
i.e. ORPHEUS conventions (a) and (b) coincide under B&G's within-group
separability ansatz, and B&G's spelling is (a)'s shape carrying (b)'s
values (Source E.1). Hébert's plain-average adjoint carrier (3.118) is
the flat-flux approximation of B&G's flux-weighted-average carrier.

---

## Missing sources — acquisition list (extraction STOPPED per brief)

Not in `scratch/literature/`. Since the 2026-07-26 B&G acquisition
these are CORROBORATION-tier (independent second sources for the E.4
χ-factorization synthesis and modern restatements), no longer blocking:

1. ~~Bell & Glasstone Ch. 6~~ — **ACQUIRED + EXTRACTED (Source E)**.
2. **Stacey, *Nuclear Reactor Physics* (2nd ed., Wiley 2007), Ch. 13**
   "Perturbation and Variational Methods" (the local Ch. 9's own adjoint
   pointer targets it) + the multigroup-constants chapter for the
   flux-weighted baseline.
3. **M.L. Williams, "Perturbation Theory for Nuclear Reactor Analysis",
   in Y. Ronen (ed.), *CRC Handbook of Nuclear Reactors Calculations*,
   CRC Press (1986)** — commonly cited from Vol. III [handbook existence
   + volume verified: OSTI 5707826; the definitive chapter-length
   treatment incl. depletion/spectral-fine-structure perturbations].
4. **Duderstadt & Hamilton, *Nuclear Reactor Analysis*, Wiley (1976)** —
   few-group-constants chapter (flux-weighted baseline + its stated
   limitations).
5. **A.F. Henry, *Nuclear-Reactor Analysis*, MIT Press (1975)** — the
   textbook home of the Henry kinetics weighting and group-collapse
   consistency.
6. The **Dorning-§8.5.4 primary chain** for the asymptotic bilinear
   result (all DOIs verified above): Larsen 1975 (JMP 16:1421), Larsen
   1976 (NSE 60:357), Zhang–Rizwan-uddin–Dorning 1995 (NSE 121:226) +
   1997 Part II (TTSP 26:433) + 1997 multiple-scales (TTSP 26:765).
   The NSE items are in the user's complete NSE archive.
7. (If the two-flux scattering rule needs a modern worked statement:
   the generalized-condensation line — Rahnema et al. — is on record in
   memory `energy-condensation-collapse-formulas` as flux-weighted GEC,
   NOT adjoint-weighted; an adjoint-weighted-condensation NSE/ANE paper
   would need a fresh targeted search once the user approves Tier-2
   extraction.)

**Zotero**: MCP tools not exposed this session — the user's annotations
could not be checked. If Bell & Glasstone / Stacey / Williams are in
Zotero with highlights on the collapse equations, those annotations
should adjudicate notation before ORPHEUS commits (AGENT.md §5
notation-oracle discipline).

---

*Memo complete. All local equations spot-verified against rendered
pages (Hébert PDF pp. 11, 19–21; Century PDF pp. 200–201, 398–399,
443–444; B&G PDF pp. 290, 297, 323–326). Verbatim quotes kept ≤3 lines
and page-cited; everything else paraphrased. Synthesis steps (E.1
carrier arithmetic, E.4 χ-factorization, E.5 rate-trade-off) are
explicitly flagged as derived-from-their-equations, not their prose.*
