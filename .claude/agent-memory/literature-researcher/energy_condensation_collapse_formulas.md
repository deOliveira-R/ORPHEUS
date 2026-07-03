---
name: energy-condensation-collapse-formulas
description: Authoritative energy-group condensation/collapse formulas (vectors, 2-axis scattering, chi) with exact eq numbers from Hebert Ch3 + Stammler Ch6, the nested-vs-non-nested verdict (MALOCS requires nesting; GROUPR/OpenMC re-integrate the continuum), the NJOY IWT within-group flux taxonomy, and the GEC least-squares-projection precedent. For ORPHEUS P5 #274.
metadata:
  type: project
---

# Energy condensation onto (non-)coincident group structures

Full memo: `.claude/plans/p5_condensation_literature.md`. ORPHEUS P5
(#274) collapses the 421-group library onto WIMS-69/172 (which are
**non-nested** — 19 boundary mismatches). Code uses `PetrovGalerkinFrame`
+ one-hot `IndicatorBasis` (= the NESTED/MALOCS case only).

## The collapse formulas — confirmed by TWO independent textbooks [both LOCAL]

**Hebert (2009) Applied Reactor Physics §3.5, pp.82-87** (NOT §4 — the
brief said §4; condensation is §3.5. Ch.4 = resonance, not local):
- (3.96) avg of a DISTRIBUTION = ∫_g du X(u) (a rate, no 1/Δu); (3.97)
  avg of a FUNCTION (flux) = (1/Δu)∫du X (lethargy avg). The key
  distinction.
- (3.103) vector: Σ_g = (1/φ_g)⟨Σφ⟩_g = (∫_g Σφ dE)/(∫_g φ dE) — flux-wt.
- (3.104) scattering 2-axis: Σ_{s,ℓ,g←h} = (1/φ_h)⟨Σ_{s,ℓ}φ⟩_{g←h};
  numerator (3.101) = ∫_g du ∫_h du' Σ(u←u')φ(u') → SINK g summed,
  SOURCE h flux-averaged, denom = φ_h (the SOURCE group). 
- (3.105) νΣ_f,g = (1/φ_g)⟨νΣ_f φ⟩_g.
- (3.112) χ_{j,g} = ∫_g du χ_j(u) — PLAIN birth-group integral, NO flux
  weight (χ integrates to 1; flux-wt would break Σχ=1). Multi-fissile
  isotope combination is (3.79)/(3.80), a DIFFERENT axis.
- Total-XS subtlety: no scalar Σ_t preserves both collision rate AND
  leakage → transport correction (3.92 Σ̄=Σ−ΔΣ_tr, 3.106 MRA vs OEWA).

**Stammler-Abbate (1983) Ch.VI §1, p.193** [`Stammler(1983)Chapter6.pdf`]
— cross-verifies Hebert exactly. ATTRIBUTION: NOT Ch.IV (=collision
probabilities; Ch.IV p.106 explicitly defers "Σ_g discussed in Ch.V.5
and VI.1"):
- VI(6a) φ_g=∫_g φ dE, χ_g=∫_g χ dE (plain). VI(6b) νΣ_f,g flux-wt.
- VI(6c) Σ_{s0,g'→g} = (∫_{g'}∫_g Σ_{s0}φ(E')dE'dE)/φ_{g'} — denom=SOURCE
  φ_{g'} (≡ Hebert 3.104). VI(6d) P1 channel is CURRENT-weighted /J_{g'}
  (a subtlety Hebert folds into the transport correction).

## Q2 nesting verdict — split by STAGE (the headline)

- **Multigroup→fewer-group (AMPX/MALOCS; ORPHEUS's stage): REQUIRES
  NESTING.** MALOCS input = a fine→broad CORRESPONDENCE ARRAY (manual
  notation `4r1 4r2 4r3` = 4 fine groups→broad 1, etc). Each fine group
  assigned ENTIRELY to one broad group ⟹ coarse boundaries MUST be a
  subset of fine boundaries; a straddling fine group is inexpressible.
  [SCALE 6.3.x manual §11.6 MALOCS2.html (host blocks fetch; via search);
  OSTI 3002301; corrob. Ann.Nucl.E 2019 DOI 10.1016/j.anucene.2019.06.025
  OSTI 1437912].
- **Pointwise→multigroup (NJOY GROUPR, OpenMC mgxs, MC²-3, Serpent): NO
  nesting** — they integrate the CONTINUOUS σ(E)φ(E) over a union grid →
  ANY structure. GROUPR Eq.70 (NJOY manual): σ_g=(∫_g F(E)σ(E)φ(E)dE)/
  (∫_g φ dE), F=feed function (=1 vector; =ℓ-th Legendre normalized
  scatter-into-g' prob for matrix). OpenMC: same ratio, arbitrary
  structure by direct CE tally.

**Fractional-overlap re-bin** (the correct construction for non-nested
collapse, ORPHEUS must reconstruct it — no named library code does it as
a primary path): σ_G=(Σ_g f_{g,G}φ_g σ_g)/(Σ_g f_{g,G}φ_g),
f_{g,G}=(∫_{g∩G} w dE)/(∫_g w dE), Σ_G f_{g,G}=1. Nested ⟹ f∈{0,1} ⟹
reduces to (3.103).

## Q3 within-group flux w(E) — NJOY IWT taxonomy [NJOY2016 manual groupx.tex, GitHub]

IWT: 1=read-in tab; 2=CONSTANT(flat-E); 3=1/E(flat-lethargy,
slowing-down default); 4=Maxwellian+1/E+fission (CANONICAL reactor;
thermal 0.025eV join 1/E@0.1eV, fission T=1.40MeV join@820.3keV);
5=EPRI-CELL PWR; 6/8 fast-reactor; 9-12 CLAW/VITAMIN-E; −n=compute flux
on the fly; 0=resonance flux from prior pass. **Standard for
CONDENSATION (≠ library gen): the computed fine-group flux φ_g
between-groups; 1/E within a split group (IWT=4 as rigorous upgrade).**

## Q4 projection precedent — YES, established

**Generalized Energy Condensation: Rahnema-Douglass-Forget 2008 NSE
160:41, DOI 10.13182/NSE160-41** [OpenAlex-confirmed, GA Tech]. Expands
within-coarse-group flux in ORTHOGONAL functions; **zeroth moment =
the standard flux-weighted average EXACTLY**, higher moments = spectral
detail. So flux-wt avg = rank-0 piecewise-constant truncation of the
projection (exactly the user's hypothesis). Consistent-GEC =
Douglass-Rahnema 2011 ANE DOI 10.1016/j.anucene.2011.09.001. PGD-in-
energy = JCP 433(2021)110744. SPH/superhomogenization = an EQUIVALENCE
correction (post-mult factors), a different lever than the weight.

ORPHEUS's `PetrovGalerkinFrame` condensation IS this projection at rank
0: trial=`IndicatorBasis` (piecewise-const), test=`WeightedIndicatorBasis
(indicator,φ)`, project=G⁻¹Mf reproduces (Σφ_gσ_g)/(Σφ_g). Enriching the
trial basis (Legendre-in-lethargy/coarse-group) ⟹ GEC, no arch change.

## ORPHEUS recommendation (Q5)

Current one-hot `IndicatorBasis` (searchsorted) = NESTED only → on
non-nested WIMS it mis-apportions straddling fine groups. Fix =
`FractionalOverlapBasis` (.evaluate returns f_{g,G} fractions, row-sum 1;
degenerates to one-hot when nested → regression-safe). Within-group w(E)
= pluggable strategy, DEFAULT 1/E (FLAT_LETHARGY), upgrades FLAT_ENERGY /
LIBRARY_SPECTRUM(IWT4) / COMPUTED_FLUX. Keep scattering = sink-summed/
source-avg (NOT both-axes-projected = the homogenize bug), χ = pure sum
χ@T. GEC rank>0 = future issue, not now.

**Provenance**: every web cite resolved to real DOI/OSTI/manual. Lewis &
Miller (brief invited it) = NOT local, eq number NOT fabricated — Hebert
+Stammler already cross-verify. Zotero not queried (local+web sufficed).
