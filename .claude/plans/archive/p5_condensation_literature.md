# Condensation onto non-coincident group structures — authoritative method

Literature pull for ORPHEUS P5 energy condensation (#274). Answers the
five mission questions with equation numbers and the exact convention.
Every source tagged **[LOCAL]** (in `scratch/literature/`) or **[web]**.

**Headline finding (the answer the mission hinges on):** the two
production stages of group-averaging are mathematically the SAME
flux-weighted reaction-rate-preserving average, but they differ
decisively in *whether they can straddle a boundary*:

- **Pointwise → multigroup** (NJOY/GROUPR, OpenMC mgxs, MC²-3, Serpent):
  integrates the *continuous-energy* `σ(E)φ(E)` directly over **any**
  group boundaries. **No nesting constraint** — the actual cross-section
  shape inside each group is the truth, so an arbitrary structure is
  trivially supported.
- **Multigroup → fewer-group** (AMPX/**MALOCS**, the case ORPHEUS is in):
  collapses *already-discretized* fine groups via a **fine→coarse
  correspondence array**. This **REQUIRES the coarse boundaries to be a
  subset of the fine boundaries (nesting)** — a fine group cannot be
  split because the input assigns each fine group entirely to one coarse
  group. Straddling a boundary requires an *added* assumption: a
  sub-fine-group flux shape `w(E)` to apportion the straddling group.

ORPHEUS is in the second stage (collapsing the 421-group library), and
its WIMS-69 / WIMS-172 structures are **non-nested** (the draft's
`boundary_mismatch_report()` already flags 19 non-coincident
boundaries). So the mission's "hard case" is real and not covered by the
nested MALOCS formula. Recommendation in §6.

---

## Q1 — The general condensation/collapse formula (CONFIRMED, two independent textbooks)

The flux-weighted, reaction-rate-preserving average is correct and is
stated identically (different equation numbers) by both authoritative
textbooks in the local folder. The defining principle in both: **the
multigroup cross section is defined so as to preserve the reaction
rate** within the group.

### Hébert, *Applied Reactor Physics* (2009), Ch. 3 §3.5 [LOCAL: `Hebert(2009)Chapter3.pdf`]

(The user cited "§4"; the condensation formalism is actually in **§3.5
"Multigroup discretization"**, pp. 82–87. Chapter 4 is resonance
self-shielding — NOT in the local folder, but it is *not* where the
collapse formulas live.)

The averaging distinguishes a **distribution** (a cross-section-weighted
quantity = a reaction rate) from a **function** (the flux):

- **Eq. (3.96)** — average of a *distribution*:
  `⟨X⟩_g = ∫_{u_{g-1}}^{u_g} du X(u) = ∫_{E_g}^{E_{g-1}} dE X(E)`
  (a plain integral — it is a rate, NOT normalized by Δu).
- **Eq. (3.97)** — average of a *function* (the flux):
  `⟨X⟩_g = (1/Δu_g) ∫ du X(u) = (1/ln(E_{g-1}/E_g)) ∫ (dE/E) X(E)`
  (the lethargy average).
- **Eq. (3.100)** — the reaction rate: `⟨Σ(r)φ(r)⟩_g = ∫_g du Σ(r,u)φ(r,u)`.

The collapse formulas (defined to preserve reaction rates), p. 84:

- **Eq. (3.103)** — vector cross section (total, etc.):
  `Σ_g(r) = (1/φ_g(r)) ⟨Σ(r)φ(r)⟩_g`
  i.e. `Σ_g = (∫_g Σ(E)φ(E)dE) / (∫_g φ(E)dE)` — **the flux-weighted average**.
- **Eq. (3.105)** — fission production: `νΣ_{f,j,g}(r) = (1/φ_g(r)) ⟨νΣ_{f,j}(r)φ(r)⟩_g`.

**Scattering matrix — the 2-axis collapse** (this is the load-bearing
one for the mission):

- **Eq. (3.101)** — the 2-axis numerator:
  `⟨Σ_{s,ℓ}(r)φ(r)⟩_{g←h} = ∫_{u_{g-1}}^{u_g} du ∫_{u_{h-1}}^{u_h} du' Σ_{s,ℓ}(r,u←u')φ(r,u')`
  (sink group `g` integrated, source group `h` integrated **with the
  flux weight on the source `u'`**).
- **Eq. (3.104)** — the collapsed matrix element:
  `Σ_{s,ℓ,g←h}(r) = (1/φ_h(r)) ⟨Σ_{s,ℓ}(r)φ(r)⟩_{g←h}`
  → **normalized by the SOURCE-group flux `φ_h`**. This is exactly the
  user's "source group flux-averaged, sink group summed": the sink axis
  `g` is summed over (the `∫_g du` runs free, accumulating into the one
  coarse sink), the source axis `h` is flux-averaged (numerator weighted
  by `φ(u')`, denominator `φ_h`).

**χ fission-spectrum collapse:**

- **Eq. (3.112)** — `χ_{j,g} ≡ ⟨χ_j⟩_g = ∫_{u_{g-1}}^{u_g} du χ_j(u)`
  → a **plain birth-group integral, NO flux weight**. χ is a probability
  distribution that integrates to 1 over all energy (Eq. 3.56), so
  collapsing it is just a sum of the probability mass landing in coarse
  group `g`. (Flux-weighting χ would destroy `Σ_g χ_g = 1`.)
  - Caveat: when collapsing a *production-weighted* χ for a multi-fissile
    mixture (combining isotopes), the average fission spectrum is the
    fission-rate-weighted combination — **Eq. (3.79)** (prompt) /
    **(3.80)** (delayed) / **(3.75)–(3.78)** define this *isotope*
    combination. But the *energy-group* collapse of an already-combined
    χ is the plain integral (3.112). ORPHEUS already separates these:
    `production-weighted chi_mix` (commit `50aa862`) handles the isotope
    axis; energy collapse is the (3.112) sum.

### Stamm'ler & Abbate, *Methods of Steady-State Reactor Physics in Simplified Geometry* (1983) [LOCAL: `Stammler(1983)Chapter6.pdf`]

**Important attribution correction:** the user cited Stamm'ler for the
collapse formula. The group-averaging definition is **NOT in Chapter IV**
(`Stammler(1983)Chapter4.pdf` — that chapter is "Integral Transport
Theory; Collision Probabilities", confirmed by reading p. 105). Ch. IV
p. 106 §2 explicitly defers: *"the definition of the group-averaged
total XS, Σ_g, is discussed in **Chapters V.5 and VI.1**."* The local
folder has **Chapter VI**; the formulas are on **p. 193 (Ch. VI §1)**:

- **Eq. VI(6a)** — energy-integrated flux/current/χ:
  `φ_{mg}(r) = ∫_g φ_m dE`; `φ_g = ∫_g φ dE`; `J_g = ∫_g J dE`;
  `χ_g = ∫_g χ dE` (χ a plain integral — matches Hébert 3.112).
- **Eq. VI(6b)** — fission production (flux-weighted):
  `νΣ_{f,g} = (∫_g νΣ_f(E')φ(E')dE') / φ_g`.
- **Eq. VI(6c)** — P0 scattering matrix (2-axis, **source-group
  flux-weighted, sink summed**):
  `Σ_{s0,g'→g} = (∫_{g'} ∫_g Σ_{s0}(E'→E)φ(E')dE'dE) / φ_{g'}`
  (note the denominator is `φ_{g'}` = the SOURCE group; identical
  structure to Hébert 3.104).
- **Eq. VI(6d)** — P1 (anisotropic) scattering component is
  **current-weighted, not flux-weighted**:
  `Σ_{s1,g'→g} = (∫_{g'} ∫_g Σ_{s1}(E'→E)J(E')dE'dE) / J_{g'}`.

**The subtlety both texts flag** (p. 193 bottom, and Hébert §3.4):
a rigorous **total** cross section is weighted by the *angular* flux and
is therefore direction-dependent — there is no single scalar `Σ_t` that
preserves both the collision rate and the leakage. Stamm'ler calls
Eq. (6d) a "nonsensical vector division" (it only makes sense if all
components of `J` share one energy shape). Hébert resolves it with the
**transport correction** (Eq. 3.92 `Σ̄ = Σ − ΔΣ_tr`, Eq. 3.106 the
condensed transport correction with the **MRA** vs **OEWA** choice).
→ For ORPHEUS: the scalar-`Σ_t` flux-weighted collapse (3.103) is the
standard pragmatic choice; the P1 anisotropic channel should ideally be
current-weighted (6d), but flux-weighting it is the common simplification.

### Lewis & Miller, *Computational Methods of Neutron Transport* (1984)

**NOT in the local folder.** Their group-collapse treatment (the same
flux-weighted reaction-rate-preserving definition, in their multigroup
chapter) would be a third corroborating textbook, but I did not verify
it against a copy — **do NOT cite a specific Lewis & Miller equation
number without acquiring it.** The two independent textbooks above
(Hébert + Stamm'ler) already cross-verify the formula exactly, so the
Q1 answer is solid without it. If the user wants the L&M citation:
"not in local folder; acquire it, or will you add it?"

### Cross-verification verdict (Q1)

Hébert (3.103)/(3.104)/(3.105)/(3.112) ≡ Stamm'ler VI(6a–6d), modulo
notation. **The general collapse formula is confirmed by two independent
authoritative sources.** Both define it to preserve the reaction rate;
both flux-weight vectors and the scattering *source* axis, sum the
scattering *sink* axis, and integrate χ without flux weight; both note
the angular-flux subtlety for the total cross section.

---

## Q2 — The non-coincident / non-nested case (THE core question)

### Do production codes require nesting? — split by stage

**Codes that collapse already-multigroup data REQUIRE nesting:**

- **AMPX / MALOCS** ("Miniature AMPX Library Of Cross Sections", the
  SCALE module; MALOCS2 in SCALE ≥6.2) [web: SCALE 6.3.x manual
  §11.6 `MALOCS2.html`; OSTI SCALE 6.3.3 manual OSTI 3002301 /
  DOI 10.2172/3002301; RSICC PSR-315]. MALOCS collapses fine-group AMPX
  master libraries into broad-group master libraries. **The broad-group
  structure is specified by a fine→broad *correspondence array*** — e.g.
  the manual's example notation `4r1 4r2 4r3` means "the first 4 fine
  groups → broad group 1, the next 4 → broad group 2, …" (the `Nr`
  repeat-count notation reads "N fine groups repeat into broad group `r`").
  Because each fine group is assigned **entirely to one broad group**,
  **MALOCS structurally REQUIRES the broad-group boundaries to be a
  subset of the fine-group boundaries (strict nesting).** A fine group
  cannot straddle a broad boundary — the input format cannot express it.
  MALOCS2 reads the collapsing spectrum from an XSDRNPM output flux file
  (or uses a built-in weighting) and flux-weights the collapse.

  > Provenance note: the `MALOCS2.html` and `t2.lanl.gov` pages are
  > indexed and real but **hard-block automated fetch** (404/empty to
  > WebFetch + curl, even with a browser UA). The `4r1 4r2 4r3`
  > correspondence-array fact is from the SCALE-manual MALOCS2 example as
  > returned by web search of that page. The structural conclusion
  > (nesting required) follows necessarily from the correspondence-array
  > input format and is independently corroborated by the AMPX/SCALE
  > fast-reactor processing paper (Jang/Kim/Lee et al.,
  > *Ann. Nucl. Energy* 2019, DOI 10.1016/j.anucene.2019.06.025 —
  > "The SCALE/AMPX multigroup cross section processing for fast reactor
  > analysis", OSTI 1437912) which describes generating a *fine* group
  > structure precisely so that broad structures collapse cleanly from it.

**Codes that go pointwise → multigroup do NOT require nesting (they
re-integrate the continuum):**

- **NJOY / GROUPR** [web: NJOY2016 manual, `github.com/njoy/NJOY2016-manual`
  `groupx.tex`; the manual report is LA-UR-17-20093, Macfarlane,
  Muir, Boicourt, Kahler, Conlin 2017]. GROUPR computes each group's
  cross section by integrating the *continuous-energy* pointwise data
  (from RECONR/BROADR) over a **union grid** of σ(E), φ(E), and the feed
  function — so it produces **any** arbitrary group structure directly.
  The manual: *"it is always possible to read in an arbitrary group
  structure … The energies must be given in increasing order."* No
  nesting concept exists because GROUPR never re-bins grouped data; it
  always re-integrates the continuum. (GROUPR's "first condensation" is
  exactly the library-generation step Hébert §3.5 references as
  occurring "in the GROUPR module, during the NJOY processing step".)

- **OpenMC `mgxs` module** [web: OpenMC docs, MGXS Part I] — computes
  multigroup XS for an **arbitrary** energy structure directly by
  tallying continuous-energy reaction rates and flux in the Monte Carlo
  run. Definition used (their Eq.):
  `σ_{x,g} = (∫∫ σ_x(r,E)Φ(r,E)dE dr) / (∫∫ Φ(r,E)dE dr)` — identical to
  Hébert (3.103). *"This approach allows arbitrary energy group
  structures without dependence on predefined fine-group libraries."*
  Same for MC²-3 (ANL, ultra-fine 2082-group → arbitrary fast structures)
  and Serpent.

**Why the asymmetry:** the pointwise codes have the *truth* (the actual
`σ(E)φ(E)` shape), so any boundary is integrable. The collapse codes
have only the fine-group *averages* — to split a fine group they must
*invent* a sub-group shape. That invented shape is Q3.

### The fractional-overlap re-binning formula

For codes/operations that *do* handle arbitrary (non-nested) structures
by re-binning already-grouped data, the apportioned collapse is:

```
            Σ_g f_{g,G} φ_g σ_g
   σ_G  =  ─────────────────────
              Σ_g f_{g,G} φ_g
```

where `f_{g,G} ∈ [0,1]` is the **fraction of fine group `g`'s reaction
rate (or flux, depending on the weight model — see Q3) that falls inside
coarse group `G``**. For a *nested* structure `f_{g,G} ∈ {0,1}` and this
reduces exactly to Hébert (3.103). For a straddling fine group,
`f_{g,G}` is a genuine fraction and `Σ_G f_{g,G} = 1` (the fine group's
contribution is partitioned, not duplicated).

`f_{g,G}` is **not** simply the energy-width overlap fraction unless the
within-group flux is flat in energy — it is the **flux-weighted overlap**:

```
              ∫_{g ∩ G} w(E) dE
   f_{g,G}  =  ─────────────────
                ∫_g w(E) dE
```

with `w(E)` the assumed sub-fine-group flux shape (Q3). This is the
"smart" part: it is exactly a re-discretization of the continuous flux
*model* restricted to fine group `g`.

> **No production *deterministic-library* code I found re-bins
> already-grouped data via this fractional-overlap formula as a primary
> path** — they either (a) require nesting (MALOCS) or (b) re-integrate
> the continuum (GROUPR/OpenMC). The fractional-overlap formula is the
> *correct* construction when you are forced to map between two
> independently-defined grouped structures (ORPHEUS's situation), and it
> is the standard formula in spectral/library *interpolation* utilities
> and in ENDF→multigroup rebinning helper tools, but be aware you are
> reconstructing rather than copying a named code's algorithm here. The
> within-group flux model (Q3) is what makes it well-defined.

---

## Q3 — The within-fine-group flux model (the NJOY IWT taxonomy)

Apportioning a straddling fine group needs `w(E)`. The reference
taxonomy is the NJOY/GROUPR **weighting-function** options, indexed by
`IWT` [web: NJOY2016 manual `groupx.tex`, §"Generalized Group
Integrals" and the IWT description list]:

| IWT | Weighting flux `w(E)` |
| --- | --- |
| 1 | read-in smooth (user-tabulated) weight function |
| 2 | **constant** (flat in energy) — used for very fine dosimetry structures |
| 3 | **1/E** (flat in lethargy) — the slowing-down default; traditional for resonance integrals |
| **4** | **thermal Maxwellian + 1/E + fission spectrum** (the canonical reactor weight; user-settable join energies/temps). The figure example: thermal `kT`=0.025 eV joined to 1/E at 0.1 eV, fission `T`=1.40 MeV joined to 1/E at 820.3 keV |
| 5 | EPRI-CELL mid-life PWR spectrum (+ fusion peak) |
| 6 | thermal — 1/E — (fission+fusion), continuous breakpoints, 300 K thermal |
| 7 | reserved |
| 8 | thermal — 1/E — fast-reactor — (fission+fusion); the fast-reactor weight |
| 9 | CLAW weight function (thermal+1/E+fission+fusion, LANL 30-group) |
| 10 | CLAW with T-dependent (recalculated Maxwellian) thermal part |
| 11 | VITAMIN-E weight (ORNL-5505): .0253-eV Maxwellian < 0.414 eV, 1/E to 2.12 MeV, 1.415-MeV fission to 10 MeV, 1/E to 12.52 MeV, fusion peak to 15.68 MeV, 1/E above |
| 12 | VITAMIN-E with T-dependent thermal part |
| −n | **compute the flux on the fly with weight function n** (the self-consistent option — a fine-group infinite-medium slowing-down flux) |
| 0 | read in a resonance flux from a prior pass (`ninwt`) |

The unifying GROUPR averaging formula (the **generalized group
integral**, NJOY manual Eq. 70):

```
            ∫_g F(E) σ(E) φ(E) dE
   σ_g  =  ──────────────────────
              ∫_g φ(E) dE
```

where `φ(E)` is the weighting flux (one of the IWT options) and `F(E)`
is the **feed function**: `F=1` for a cross-section vector; `F = ℓ-th
Legendre component of the normalized probability of scattering into
secondary group g'` for the scattering matrix; `F = photon yield` for
production. (This *is* Hébert 3.103/3.104 with the feed function
selecting the channel — same math, NJOY's unifying packaging.)

### Which model is standard, and for which task

- **Library generation (pointwise → fine multigroup):** the **IWT=4 /
  IWT=−n** family (fission + 1/E + Maxwellian, or a computed
  infinite-medium flux) is standard — it is the same `w(E)` used to
  generate the fine `σ_g` in the first place. This is the NJOY GROUPR
  step Hébert §3.5 names.
- **Condensation (fine multigroup → fewer-group), the ORPHEUS case:**
  the **right weight is the actual computed fine-group flux `φ_g` from a
  transport solve** (a lattice/cell calculation). This is what MALOCS2
  does (reads the XSDRNPM flux), what Stamm'ler VI(6) assumes (`φ` is the
  solved flux), and what Hébert §3.5 means by the "next energy
  condensation in the lattice code". The point of condensation is to bake
  the *problem-specific* spectrum into the few-group constants; using a
  generic library spectrum would defeat it.
  - **But** for the *sub-fine-group* shape needed only to apportion a
    *straddling* fine group (the `f_{g,G}` fractions), the computed flux
    gives you only the group *average* `φ_g`, not the within-group shape.
    There you fall back to a within-group model: **flat-in-lethargy
    (1/E)** is the standard slowing-down default for the resonance/fast
    range; **flat-in-energy** for very narrow groups; and the **library
    weighting spectrum** (IWT=4) is the most defensible because it is
    self-consistent with how the fine `σ_g` were generated. The hierarchy
    of fidelity: computed `φ_g` (between groups) × library/1-E shape
    (within a split group).

**Standard choice for condensation, stated plainly:** flux-weight with
the **computed fine-group flux** `φ_g` between groups; for the
within-group sub-shape needed to split a straddling group, use the
**1/E (flat-in-lethargy)** slowing-down model in the epithermal/fast
range as the default, with the **library weighting spectrum (IWT=4
fission+1/E+Maxwellian)** as the more rigorous (and self-consistent)
upgrade. Flat-in-energy is acceptable only when groups are narrow enough
that the choice is immaterial.

---

## Q4 — Is condensation ever posed as a least-squares / projection? (YES — strong precedent)

There is a clear, named literature line that formulates energy
condensation as a **spectral projection onto an orthogonal-function
basis**, generalizing the simple flux-weighted average. This is directly
relevant: it *is* the user's "treat the spectrum as a measure and
project the XS onto a coarse basis" idea, and it is exactly how ORPHEUS
already implements condensation (the `PetrovGalerkinFrame`).

### Generalized Energy Condensation (GEC) — the primary precedent

- **Rahnema, Douglass & Forget (2008), "Generalized Energy Condensation
  Theory", *Nuclear Science and Engineering* 160:41–58, DOI
  10.13182/NSE160-41** [web; OpenAlex confirmed: NSE, 2008, Georgia
  Tech]. Abstract (verbatim excerpt): *"A generalization of multigroup
  energy condensation theory has been developed … This is accomplished
  by **expanding the energy dependence of the angular flux in a set of
  general orthogonal functions**. The expansion leads to a set of
  equations for the angular flux moments in the few-group framework. The
  zeroth moment generates the standard few-group equation while the
  higher-moment equations generate the detailed spectral resolution
  within the few-group structure."*
- **Douglass & Rahnema (2011), "Consistent generalized energy
  condensation theory", *Annals of Nuclear Energy*, DOI
  10.1016/j.anucene.2011.09.001** [web] — the consistent (flux-moment
  recoupling) refinement.

**Relationship to the flux-weighted average (the exact answer to "do
they coincide for piecewise-constant bases"):** YES. GEC expands the
within-coarse-group flux as `φ(E) ≈ Σ_n φ_{n,G} P_n(E)` in orthogonal
functions `P_n` over coarse group `G`. The **zeroth moment** (n=0,
i.e. the constant/piecewise-constant basis function on `G`) **recovers
the standard flux-weighted multigroup average exactly** — the paper
states this explicitly ("the zeroth moment generates the standard
few-group equation"). The higher moments (n≥1) add the within-coarse-
group spectral detail that the simple average throws away. So the
flux-weighted average is the **rank-0 / piecewise-constant truncation**
of the projection — precisely the relationship the user hypothesized.

### Other projection / reduced-order precedents

- **Proper Generalized Decomposition (PGD) in energy** — "Reduced-order
  modeling of neutron transport separated in energy by Proper
  Generalized Decomposition", *J. Comp. Phys.* 433 (2021) 110744
  [web, OpenAlex]. A low-rank separated-representation of the energy
  dependence — the modern reduced-order cousin of GEC.
- **SPH / superhomogenization equivalence factors** (Hébert §4 / Kavenoky;
  Sanchez) — not strictly a least-squares projection, but an *equivalence*
  formulation: it does not change the weighting; it post-multiplies the
  collapsed XS by factors `μ_G` that force the few-group solution to
  reproduce the reference reaction rates (a *correction* to the projection,
  used when the simple flux-weighted collapse loses reaction-rate
  consistency across the coarse interfaces). Cite if the user wants
  equivalence on top of condensation; it is a different lever than the
  weight choice.

### Verdict (Q4)

The least-squares/projection framing is **not exotic — it is an
established generalization (GEC, 2008–2011) and it degenerates to the
flux-weighted average at zeroth order on a piecewise-constant basis.**
ORPHEUS's existing `PetrovGalerkinFrame` condensation IS this projection
at rank 0: trial basis = `IndicatorBasis` (piecewise-constant on coarse
groups), test basis = `WeightedIndicatorBasis(indicator, φ)` (the
flux-weighted indicator), and `project = G⁻¹ M f` reproduces
`Σ_G = (Σ φ_g Σ_g)/(Σ φ_g)`. The clean generalization path is to enrich
the trial basis beyond rank-0 indicators (Legendre-in-lethargy per
coarse group) to get GEC — the Frame machinery already supports it.

---

## Q5 (the deliverable) — Concrete recommendation for ORPHEUS

**Context (verified against the code this session via the P5.0 gate
plan `.claude/plans/archive/p5_0_condensation_gate.md`):** ORPHEUS condensation
is built on the `PetrovGalerkinFrame` (`orpheus/numerics/frame.py:295`),
with the partition as a one-hot `IndicatorBasis`
(`indicator_basis.py`, `searchsorted` one-hot — **each fine group → one
coarse group**). This is the **nested (MALOCS-equivalent) case**. The
WIMS-69/172 structures are **non-nested** (`boundary_mismatch_report()`
flags 19 boundaries; e.g. 172-groups 1–5 sit above the 69-group
ceiling). So the current one-hot basis **cannot represent a fine group
that straddles a coarse boundary** — it will assign each straddling fine
group wholly to one side (the `searchsorted` containing-interval rule),
silently mis-apportioning its rate.

### Recommendation

1. **Use the flux-weighted reaction-rate-preserving collapse — confirmed
   canonical** (Hébert 3.103/3.104/3.105/3.112 ≡ Stamm'ler VI(6a–6d)).
   ORPHEUS already does this correctly for the nested channels: vectors
   `σ_G = project(σ)`, scattering 2-axis = **sink-summed, source-flux-
   averaged** (NOT both-axes-projected — that is the `homogenize`
   behavior and is wrong for condensation), χ = **pure birth-group sum**
   (`χ @ T`, never projected). Keep these exactly as the P5.0 gate pins
   them.

2. **For the non-nested boundary, generalize the partition basis from a
   one-hot indicator to a fractional-overlap operator.** Replace the
   `{0,1}` membership `T[g,G]` with the flux-weighted overlap fraction
   `f_{g,G}` (Q2 formula). Concretely: the partition matrix becomes
   `f_{g,G} = (∫_{g∩G} w(E)dE) / (∫_g w(E)dE)` with row-sums
   `Σ_G f_{g,G} = 1`. This is a *minimal, principled* generalization of
   the existing `IndicatorBasis` — a **`FractionalOverlapBasis`** whose
   `.evaluate` returns overlap fractions instead of a one-hot row. The
   `PetrovGalerkinFrame` machinery (`project = G⁻¹ M f`) is unchanged;
   only the basis changes. For a nested structure it degenerates to the
   current one-hot exactly (regression-safe).

3. **Within-group flux default: make it pluggable, default to 1/E
   (flat-in-lethargy).** The `f_{g,G}` integrals need `w(E)`. Provide a
   small enum/strategy:
   - `FLAT_LETHARGY` (1/E) — **the default** (slowing-down standard,
     Q3; matches NJOY IWT=3 and the universal epithermal assumption).
   - `FLAT_ENERGY` (constant) — for narrow groups (NJOY IWT=2).
   - `LIBRARY_SPECTRUM` (fission + 1/E + Maxwellian, NJOY IWT=4) — the
     rigorous self-consistent upgrade; reuses the same `w(E)` the
     421-group library was generated with.
   - (future) `COMPUTED_FLUX` — reconstruct the within-group shape from a
     finer auxiliary solve.
   Pluggable because the right choice is regime-dependent (the user's
   resonance/thermal groups where the coarse grid is *locally finer* than
   the fine grid are exactly where the within-group model matters most),
   and because it is the natural seam — the *between-group* weight stays
   the computed `φ_g`, only the *within-split-group* sub-shape is the
   plug.

4. **Document the two-stage distinction in the theory page.** ORPHEUS is
   in the *collapse* stage (fine multigroup → fewer), not the *library*
   stage (pointwise → multigroup). The nesting requirement is a property
   of the collapse stage; the fractional-overlap generalization is what
   lifts that requirement. Note that MALOCS *requires* nesting and
   ORPHEUS is deliberately going beyond it — this is a real capability
   the production deterministic-library codes mostly do not have (they
   sidestep it by re-integrating the continuum, which ORPHEUS cannot do
   from a pre-grouped 421-library).

5. **Keep the door open to GEC (Q4) as the rank>0 future.** The
   fractional-overlap rank-0 collapse is the right *now*. If within-
   coarse-group spectral fidelity is ever needed (deep-penetration,
   strong spectral gradients), the same `PetrovGalerkinFrame` enriched
   with Legendre-in-lethargy trial functions per coarse group gives
   Generalized Energy Condensation (Rahnema-Douglass-Forget 2008), of
   which the current flux-weighted average is the zeroth moment. No
   architectural change — just a richer trial basis. File as a future
   issue, do not build now.

### One correctness caveat to carry into implementation

The total cross section has **no** weight that simultaneously preserves
collision rate *and* leakage (Hébert §3.4 / Stamm'ler p. 193). The
flux-weighted `Σ_t,G` (3.103) preserves the collision rate; for the
transport/leakage-consistent value, the **transport correction** (Hébert
3.92/3.106, MRA or OEWA) is applied, and the **P1 scattering channel is
properly current-weighted** (Stamm'ler 6d), not flux-weighted. ORPHEUS's
current collapse flux-weights all channels uniformly — this is the
standard pragmatic choice, but flag in the theory page that the rigorous
P1/transport treatment differs, so a future anisotropic-condensation
upgrade knows where the approximation lives.

---

## Source ledger

**[LOCAL] (`scratch/literature/`) — read this session:**
- Hébert, *Applied Reactor Physics* (2009), Ch. 3 §3.5, pp. 82–87 —
  Eqs. (3.96)–(3.112). `Hebert(2009)Chapter3.pdf`. The primary collapse
  source. (Ch. 4 resonance/condensation NOT local; collapse formulas are
  in Ch. 3, not Ch. 4 as the brief assumed.)
- Stamm'ler & Abbate, *Methods of Steady-State Reactor Physics in
  Simplified Geometry* (1983), **Ch. VI §1, p. 193** — Eqs. VI(6a)–(6d).
  `Stammler(1983)Chapter6.pdf`. Cross-verifies Hébert. (NOT Ch. IV —
  that is collision probabilities; Ch. IV p. 106 defers to Ch. V.5/VI.1.)

**[web] — verified against a real database / authoritative manual:**
- NJOY2016 manual (Macfarlane, Muir, Boicourt, Kahler, Conlin 2017,
  LA-UR-17-20093), GROUPR chapter `groupx.tex` —
  `github.com/njoy/NJOY2016-manual` (GitHub raw, read this session).
  Eq. 70 (generalized group integral + feed function), full IWT
  weighting taxonomy.
- SCALE 6.3.x manual §11.6 MALOCS2 — `scale-manual.ornl.gov/MALOCS2.html`
  (page real & indexed; host blocks automated fetch — content via web
  search of the page). OSTI SCALE 6.3.3 manual OSTI 3002301,
  DOI 10.2172/3002301. The `Nr` correspondence-array → nesting
  requirement.
- Jang/Kim/Lee et al. (2019), "The SCALE/AMPX multigroup cross section
  processing for fast reactor analysis", *Ann. Nucl. Energy*,
  DOI 10.1016/j.anucene.2019.06.025 (OSTI 1437912) — corroborates the
  fine-structure-then-collapse workflow.
- Rahnema, Douglass & Forget (2008), "Generalized Energy Condensation
  Theory", *Nucl. Sci. Eng.* 160:41–58, **DOI 10.13182/NSE160-41**
  (OpenAlex-confirmed). The projection precedent (Q4).
- Douglass & Rahnema (2011), "Consistent generalized energy condensation
  theory", *Ann. Nucl. Energy*, DOI 10.1016/j.anucene.2011.09.001.
- "Reduced-order modeling of neutron transport separated in energy by
  PGD", *J. Comp. Phys.* 433 (2021) 110744 (OpenAlex-confirmed).
- OpenMC `mgxs` documentation (Part I) — `docs.openmc.org` —
  arbitrary-structure direct-tally definition.

**NOT verified / acquire-if-needed (do NOT cite an equation number):**
- Lewis & Miller, *Computational Methods of Neutron Transport* (1984) —
  not in local folder; their multigroup-collapse chapter would be a third
  corroboration but was not read. (Hébert+Stamm'ler already cross-verify.)
- NRC ML12184A002 (HTGR broad-group library via MALOCS) — a practical
  MALOCS application that would give a verbatim nesting statement; the
  PDF timed out on fetch. Not load-bearing (the SCALE manual + AMPX paper
  already settle nesting).

**Provenance discipline applied:** every web citation resolved to a real
DOI/OSTI-ID/manual-of-record (no memory-only or docstring-only
citations). The MALOCS nesting conclusion is from the manual's input
format + an independent corroborating paper, not a single source. The
"Lewis & Miller §X" the brief invited was deliberately NOT fabricated —
flagged as unacquired.

**Zotero:** not queried this session — the local folder + web databases
covered every question, and the two load-bearing collapse equations are
in local PDFs whose extraction needed no annotation oracle. If the user
wants their Zotero highlights on the Rahnema-Douglass-Forget GEC paper or
the SCALE manual surfaced, re-run with Zotero up.
