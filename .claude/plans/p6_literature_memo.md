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

**Status**: COMPLETE (2026-07-25). Three local sources extracted +
verified; per-channel classical formulas traced to five missing classics
(acquisition list at end); no local-source conflicts.

---

## 0. Local-library inventory verdict

Grep sweep of all sidecars for `adjoint`, `bilinear`, `perturbation`,
`condensation`, `homogeniz*`, `group constant`, `collapse`:

**FOUND LOCALLY (extractable content):**

| File | Relevant content |
|---|---|
| `Hebert(2009)Chapter3.pdf` | Adjoint transport equation, Rayleigh ratio + first-order stationarity statement, multigroup adjoint flux, flux-weighted condensation §3.5 (the baseline the bilinear collapse corrects) |
| `Nuclear Computational Science - A Century in Review.pdf` | Ch. 4 "Reactor Core Methods" §4.3.2 homogenization: importance-weighted reaction-rate conservation ⟨Φ†,Q⟩, flux-volume average, correction factors, GET/discontinuity factors; a later chapter (§8.5.4) explicitly states the asymptotic-homogenization result that homogenized constants are weighted by BOTH forward and adjoint flux |

**NOT IN LOCAL FOLDER (the classical per-channel bilinear collapse
prescriptions live in these — extraction STOPPED on each per the brief;
acquisition list at the end):**

- Bell & Glasstone, *Nuclear Reactor Theory* (1970) — Ch. 6 (adjoint,
  perturbation theory, variational methods; the classic multigroup
  adjoint-weighted group-constants discussion §6.4).
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

## Consolidated taxonomy table

Channel class → classical prescription → where it is grounded LOCALLY →
status. "Bilinear pairing" below = the duality contraction the collapse
must preserve; the per-channel fine→coarse COLLAPSE FORMULAS themselves
are in the missing classics (last column).

| # | Channel class | Classical (bilinear / eigenvalue-consistent) prescription | Local grounding (verified) | Per-channel formula w/ eq. number |
|---|---|---|---|---|
| 1 | Vector channels (Σ_t, Σ_a, Σ_f, νΣ_f) | Preserve the adjoint-weighted reaction worth: Σ_G from ⟨φ*Σφ⟩_G/⟨φ*φ⟩_G-type ratios so that δk of the collapse vanishes at first order | Stationarity of K_eff w.r.t. δΣ at fixed φ,φ*: Hébert (3.63)+p. 77 statement. "Weighted by both the transport theory flux and the … adjoint flux" for homogenized constants: Dorning §8.5.4 p. 436 (from Fredholm solvability, refs [77]–[81]) | NOT LOCAL — B&G Ch. 6 / Stacey Ch. 13 / Williams CRC / Henry. Flux-weighted baseline (3.103) is local |
| 2 | Scattering matrix Σ_s(g′→g) | Two-flux rule: each (g′,g) pair weighted by SOURCE forward flux φ_{g′} AND SINK adjoint φ*_g (the pairing in ⟨φ*,Sφ⟩ = Σ_g Σ_{g′} φ*_g Σ_{s,g←g′} φ_{g′}) | The bilinear pairing structure is explicit in Hébert (3.61) continuous / (3.121) multigroup (transposed kernel: sink-adjoint contraction) and Dorning (8.30)+(8.37). No local source writes the fine→coarse two-flux collapse | NOT LOCAL — same classics. Flux-weighted baseline (3.104): source-flux-averaged, sink-summed, φ* absent |
| 3 | Fission source: χ and νΣ_f | Factored rank-1 treatment: adjoint contracts the EMISSION spectrum (⟨χφ*⟩), forward contracts PRODUCTION (⟨νΣ_f φ⟩); preserved invariant = total fission worth ⟨φ*,Fφ⟩ | Rank-1 factorization is explicit in Hébert (3.61), (3.63) numerator = Σ_j ⟨χ_j φ*⟩·⟨νΣ_{f,j}φ⟩, and Dorning F_t (8.39)–(8.41); worked collapsed-spectrum-channel example = β̄ᵢ (8.45) = (N₀†,M_tⁱΨ)/F_t | NOT LOCAL for the explicit χ_G formula. Flux-weighted baseline: χ_{j,g} plain sum (3.112); kinetics χ^pr/χ^del collapse (3.132)–(3.133) weighted by νΣ_fφ only |
| 4 | Perturbation-theory grounding (δk formula) | δρ = ⟨φ*, δH φ⟩/⟨φ*, Fφ⟩; bilinear-weighted constants zero the first-order collapse error because K_eff is stationary in (φ,φ*) | Dorning (8.44) ρ = (N₀†,(H_t−H₀)Ψ)/F_t with F_t = (N₀†,M_tΨ) — the formula EXACTLY, reactivity flavor; Hébert (3.63) + the p. 77 stationarity statement ("δK_eff … without using δφ or δφ*") | LOCAL (both). The formal "collapse error is second order" THEOREM statement: not local (B&G/Williams) |
| 5 | Condensation (energy) vs homogenization (space) | Energy condensation: classical texts give flux-weighted collapse + bilinear correction from perturbation theory. Spatial homogenization: EITHER flux-weight + equivalence factors (μ/SPH, discontinuity factors) OR asymptotic multiscale which FORCES forward×adjoint weighting | Roy §4.3.2 (4.41)–(4.47): industrial spatial chain = flux-volume average (4.44) + μ factors + GET f^± (4.47); adjoint only for response functionals ⟨Φ†,Q⟩. Dorning §8.5.4: asymptotic route yields bilinear weights + discontinuity factors at every scale (2- and 3-scale) | Hébert §3.5 = the flux-weighted energy-condensation baseline (LOCAL). The explicit "what flux weighting gets wrong at first order" discussion: NOT LOCAL |
| 6 | (bonus) Kinetics parameters β_eff, Λ | β̄ᵢ = ⟨φ*,M^iφ⟩/⟨φ*,Fφ⟩, Λ̄ = ⟨φ*,(1/v)φ⟩/⟨φ*,Fφ⟩ — the oldest adjoint-weighted collapse (to a 0-D model) | Dorning (8.45)/(8.46) [LOCAL, verified]; lineage Ussachoff 1955 / Henry 1955, 1958 | LOCAL in operator form. Hébert's lattice-code (3.126)–(3.135) are the flux-weighted counterparts — a live example of the two conventions coexisting |
| 7 | (bonus) Multigroup-adjoint consistency | "Transpose the multigroup operator" ≠ "discretize the continuous adjoint" unless constants are consistently weighted; adjoint carrier is a group-AVERAGE (function), forward a group-INTEGRAL (distribution) | Hébert (3.118)–(3.121) + the distribution/function statements pp. 77, 85 [LOCAL, verified]. Hébert does NOT state the consistency condition | The consistency-condition discussion: NOT LOCAL (B&G Ch. 6 is its classic home) |

**Agreement/disagreement between the local sources**: no conflicts —
the three local sources partition cleanly: Hébert = adjoint objects +
flux-weighted collapse baseline; Roy = industrial spatial homogenization
(flux weight + equivalence corrections, adjoint for responses only);
Dorning = the bilinear prescriptions (kinetics collapse exactly; spatial
homogenization as a solvability RESULT). The known literature
DISAGREEMENT the brief asks about (per-pair sink×source vs source-only
scattering weighting) cannot be adjudicated from the local corpus: the
local sources contain the bilinear PAIRING structure but no source
states the fine→coarse scattering collapse rule.

---

## Missing sources — acquisition list (extraction STOPPED per brief)

Not in `scratch/literature/`; needed to pin the per-channel classical
collapse formulas with equation numbers. In priority order:

1. **Bell & Glasstone, *Nuclear Reactor Theory*, Van Nostrand Reinhold
   (1970), Ch. 6** ("The Adjoint Equation, Perturbation Theory, and
   Variational Methods") — the classic multigroup adjoint + consistent
   group-constants discussion (the taxonomy rows 1, 2, 4, 7). The single
   highest-value acquisition.
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
443–444). Verbatim quotes kept ≤3 lines and page-cited; everything else
paraphrased.*
