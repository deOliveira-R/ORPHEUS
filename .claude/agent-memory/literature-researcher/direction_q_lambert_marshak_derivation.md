---
name: Direction Q — Lambert/Marshak primitive mismatch literature pass
description: Pursuit of principled derivation of F.4's Lambert (no µ weight on outgoing) + Marshak (µ-weighted) primitive mismatch. Sanchez 2014 NSE 177 abstract retrieved in full — gives a rigorous degeneracy theorem for first-order P_N that allows multiple admissible IC/BC families but does NOT address IC-method escape/transmission conventions. Bogado Leite 1998 ANE 25 remains closed-access and unreferenced even in its author's own later work. Corngold 2002 + 2004 addendum, Wio 1984, Krishnani 1982/1985 all behind Elsevier paywall with no abstract anywhere.
type: project
---

# Direction Q — Lambert/Marshak primitive-mismatch literature pass

**Author**: literature-researcher agent, 2026-04-22.
**Task**: GitHub Issue #122 — find a principled origin for the
Lambert (sin θ · exp(−τ), no outgoing µ weight) + Marshak
(µ-weighted W) primitive mismatch that is empirically load-bearing
in F.4.
**Budget spent**: ~55 min literature (Tier-2 heavy, Zotero
confirmed in its flaky mode).
**Verdict in one line**: Sanchez 2014 gives a gauge-theoretic
principle ("non-uniqueness of first-order P_N IC/BC") but for
*differential* P_N, not for integral-CP primitives. Bogado Leite
1998 remains a locked door — abstract invisible, no OA copy,
authoring ecosystem shows zero build-on. Corngold / Wio / Krishnani
all locked the same way. **Direction Q cannot be closed from public
metadata alone; a PDF of Bogado Leite 1998 via institutional ILL is
the single highest-value follow-up.**

---

## 1. Status per reference

| Reference | DOI | OA? | PDF access | Abstract | Grade |
|-----------|-----|-----|------------|----------|-------|
| **Sanchez 2014 NSE 177(1) 19–34** | 10.13182/NSE12-95 | No (Unpaywall confirmed `oa_status: closed`, `has_repository_copy: false`) | Not found — HAL empty, CEA archive empty, ResearchGate blocks, Tandfonline 403 | **Full abstract retrieved from OpenAlex (reproduced below §2.1)** | **Abstract only** |
| **Bogado Leite 1998 ANE 25(1–3) 129–139** | 10.1016/S0306-4549(97)00026-1 | No | Not found — author's own PPGEM thesis group does not cite it; CNEN has no OA deposit | **None from any API (CrossRef/OpenAlex/S2 all empty)** | **Not found** |
| **Bogado Leite 1999 TTSP 28(4)** | 10.1080/00411459908206038 | No | Not found | Abstract retrieved (reproduced §2.2) — this is the 1999 follow-up that reports "*significant discrepancies*" between first-flight probability formulations, but does not touch the interface-current revision | Abstract only |
| **Corngold 2002 ANE 29(5) 509–523** | 10.1016/S0306-4549(01)00102-5 | No | Not found — Caltech archive blocks web crawl | **None** | **Not found** |
| **Corngold 2004 addendum** | 10.1016/j.anucene.2003.11.005 | No | Not found | **None** | **Not found** |
| **Wio 1984 ANE 11(11) 561–?** | 10.1016/0306-4549(84)90061-6 | No | Not found | **None** | **Not found** |
| **Krishnani 1982 ANE 9(5)** | 10.1016/0306-4549(82)90083-4 | No | Not found | **None** | **Not found** |
| **Krishnani 1985 ANE 12(**?**)** | 10.1016/0306-4549(85)90080-5 | No | Not found | **None** | **Not found** |

Note on institutional PDF procurement: the only OA follow-up
associated with Bogado Leite is a UFRGS master's thesis by Letícia
Jenisch Rodrigues (2011, `hdl.handle.net/10183/28965`). I fetched
its full text. The thesis is a PIJM collision-probability and
Dancoff-factor code report; it lists Bogado Leite's 2000 / 2001 /
2003 / 2004 papers but **does not cite his 1998 ANE paper** at all.
That fact compounds the orphan status: not even his co-advisee's
PhD-adjacent work picks up the 1998 revision.

---

## 2. Key excerpts (verbatim from public APIs)

### 2.1 Sanchez 2014 — full abstract (via OpenAlex W-record for DOI 10.13182/NSE12-95)

Affiliation on file (via Unpaywall `z_authors`): *CEA de Saclay,
Laboratoire du Transport Stochastique et Déterministe, DEN/DM2S/SERMA
91191, Gif-sur-Yvette cedex, France.*

> We investigate the degeneracy of the first-order $P_N$ equations
> and construct interface and boundary conditions that ensure a
> unique solution. Our technique is based on establishing an
> equivalence between the first- and second-order $P_N$ equations
> and showing that the (regular) second-order equations with
> opposite parity to $N$ are nondegenerate. Assuming bounded
> angular flux moments and sources, we derive interface and
> boundary conditions for the regular second-order equations that,
> via the equivalence, are those to be used with the first-order
> $P_N$ equations. While providing independent derivations, our
> results reproduce those derived using solid harmonic expansions
> by Davison and Rumyantsev in the 1950s.

**What this tells us** (no speculation beyond the abstract):

1. Sanchez identifies a *degeneracy* of the first-order $P_N$
   equations. "Degeneracy" here means the first-order equations
   alone do not determine a unique solution — additional IC/BC
   closures are needed, and the choice is not uniquely constrained
   by the first-order equations themselves.

2. His technique: *build the equivalent second-order $P_N$ system
   with opposite parity*, and derive IC/BC from the second-order
   equations (which are nondegenerate). Pull back to first-order
   via equivalence.

3. His derivation *reproduces Davison and Rumyantsev (1950s)*, who
   used solid-harmonic expansions. The solid-harmonic machinery is
   the older, equivalent path.

**What this does NOT tell us** — this is the critical load-bearing
limit of the abstract: Sanchez 2014 is *about $P_N$ approximations
of the differential transport equation*, not about the *integral
CP/IC equations* in which F.4 lives. F.4's Lambert-vs-Marshak
asymmetry is a choice of *angular weight* on the **escape /
transmission response matrices** $P_\mathrm{esc}$, $G_{bc}$, $W$ of
the integral-CP method, not a choice of moment closure on the
first-order $P_N$ flux expansion. Sanchez's degeneracy theorem
applies to the left-hand side ODE system; F.4's primitives are
quadratic forms on the right-hand side integral operator.

**Relevance verdict**: **partially addressed**. Sanchez 2014 makes
a *philosophical* case that IC/BC choice is gauge-like for
first-order $P_N$ — the equations themselves do not impose a
unique closure. This is a precedent for the *idea* that ORPHEUS's
"Lambert + Marshak" combination might be an admissible gauge for
its own IC formulation. But the abstract and the derivation scope
(ODE $P_N$) do not supply the solid-harmonic identity that would
*tell us why* Lambert P/G + Marshak W cancels the rank-1 error and
fails at rank-2+. We do not have access to the main text equations.

**Single most useful fact recoverable from this paper without a
PDF**: the paper *cites Davison (1957) and Rumyantsev (mid-1950s
Soviet literature)* for an earlier solid-harmonic argument. That
earlier material is what a principled derivation would build on.
Sanchez 2014 is the modern re-expression. **The pointer to
Davison/Rumyantsev is the signal; the 2014 paper itself is a
synthesis we cannot use without the body text.**

---

### 2.2 Bogado Leite 1999 — full abstract (TTSP 28(4), DOI 10.1080/00411459908206038)

Author affiliation on record: CNEN/RJ (Brazilian National Nuclear
Energy Commission). OpenAlex institution: "National Nuclear Energy
Commission."

> **Abstract.** First-flight probabilities for lattice integral
> transport calculation, reported in the literature under different
> assumptions, are compared. Errors on calculated escape and
> collision probabilities in cells with variable number of annular
> regions are assessed, as well as differences in the Dancoff
> factor approximations, in cylindrical and spherical geometries.
> Numerical examples for selected situations show, in some cases,
> significant discrepancies among the results.

**What this tells us**:

1. By 1999 Bogado Leite is comparing *first-flight probability
   formulations across the literature* in cylindrical and spherical
   cells. He explicitly notes "*significant discrepancies among the
   results*" for the same lattice problem under different
   formulations. This is precisely the phenomenology the ORPHEUS
   team is confronting at rank-N.

2. The title of the 1998 paper ("**Revised** interface-current
   relations...") combined with the 1999 paper's tone implies that
   Bogado Leite's 1998 paper is proposing *one specific revision*
   (and hence a choice of convention) among several in
   circulation — and that the choice matters numerically.

**What this does NOT tell us**: which conventions he compared,
which he preferred, what the "revision" in the 1998 paper actually
changes relative to the baseline (Stamm'ler / Carlvik / Roy), or
whether his "revision" maps onto either the Lambert or the Marshak
side of ORPHEUS's asymmetry. For that we need the 1998 main text.

**Relevance verdict**: **partially addressed**. The 1999 abstract
is direct evidence that a *convention-choice inconsistency* in
cylindrical/spherical first-flight probabilities was a known source
of "significant discrepancies" as of 1999, and that Bogado Leite
had published a "revised" formulation one year earlier. But we do
not have the equations that would tell us whether his revision
coincides with, contradicts, or sheds light on F.4's Lambert +
Marshak pairing.

---

### 2.3 Corngold 2002 + 2004 addendum, Wio 1984, Krishnani 1982 + 1985

**No abstracts available from any public API** (OpenAlex, CrossRef,
Semantic Scholar all return empty `abstract` fields). All are
Elsevier-era ANE papers from 1982–2004 that predate the 2008-ish
abstract-deposition norm. The metadata available:

- **Corngold 2002** — N. R. Corngold, "On the collision probability
  for the infinite cylinder," *Annals of Nuclear Energy* (Pergamon /
  Elsevier), DOI 10.1016/S0306-4549(01)00102-5. Cited by 11, of
  which 1 is classified "influential" by Semantic Scholar. Title
  alone makes this paper the canonical modern derivation of $P_c$
  in infinite cylinder; standard expected content is the $Ki_3$
  Peierls integral. Abstract unavailable.
- **Corngold 2004** — Addendum to the above, DOI
  10.1016/j.anucene.2003.11.005. Zero citations (sign of a short
  corrections note). Abstract unavailable.
- **Wio 1984** — H. S. Wio, "Transformation law for response
  fluxes within the collision-probability method," ANE 11(11),
  DOI 10.1016/0306-4549(84)90061-6. Cited by 2. Abstract
  unavailable. Title-level relevance: a *response-flux
  transformation* is the closest pre-1990 title match to the
  "gauge / convention rescaling" framing we are pursuing.
- **Krishnani 1982** — P. D. Krishnani, "Interface current method
  for PHWR cluster geometry with anisotropy in the angular flux at
  interfaces," ANE 9(5), DOI 10.1016/0306-4549(82)90083-4. Cited
  by 11. Explicitly anisotropic IC for cluster geometry. Abstract
  unavailable.
- **Krishnani 1985** — P. D. Krishnani, "Transformation laws for
  probabilities in the interface current method and its application
  to cluster geometry," ANE 12(?), DOI 10.1016/0306-4549(85)90080-5.
  Cited by 2. Title-exact match for "transformation laws in IC."
  Abstract unavailable.

**Relevance verdict for each**: **not in the literature checked
(at metadata level)**. None of the five papers expose an abstract;
the 1980s deposition vintage means no abstract ever existed in
indexed form. The *titles* are all directly on point, particularly
Krishnani 1985 and Wio 1984. Without PDFs these are gray literature
to the current investigation.

---

## 3. What was verified in the author ecosystem

### 3.1 Bogado Leite — full author bibliography (OpenAlex A5038781981)

16 works total indexed under this author ID. The nuclear-relevant
subset:

| Year | Venue | Title (truncated) |
|------|-------|-------------------|
| 1998 | ANE 25 | **Revised interface-current relations for the unit-cell transport problem in cylindrical and spherical geometries** |
| 1999 | TTSP 28 | A comparison of first-flight probabilities in multiregion cylindrical and spherical cells |
| 2000 | PNE 22 | An assessment of WIMS method for computing collision probabilities in one-dimensional annular systems |
| 2000 | Kerntechnik 65 | Improvements in calculating collision probabilities in cluster geometry |
| 2001 | Kerntechnik 66 | Collision probability calculation in cluster lattices with square outer boundaries |
| 2003 | — | Determinação de Probabilidades de Colisão em Geometria Cluster (Monte Carlo) |
| 2004 | Kerntechnik 69 | Improved Dancoff Factors for Cluster Fuel Bundles by the WIMSD Code |
| 2008 | Kerntechnik 73 | Deterministic calculation of grey Dancoff factors in cluster cells with cylindrical outer boundaries |
| 2008 | ANE | Solution of the neutron transport problem with anisotropic scattering in cylindrical geometry by the decomposition method |
| 2012 | PNE | Dancoff factors with partial absorption in cluster geometry and their effects on unit cell multiplication |

**Observation**: The 1998 ANE paper's title promises a hollow-cell /
unit-cell IC revision; every subsequent Bogado Leite paper pivots to
**collision probability** (Dancoff factors, PIJM, cluster CP). No
later paper of his advertises an IC-rank-N development. This is
consistent with one of two possibilities:

(a) the 1998 paper's "revision" was a scalar-CP correction that
    needed no follow-through, and subsequent papers quietly *use*
    it under the hood in PIJM — in which case the "revision" is a
    Stamm'ler-equivalent tweak, **not** a rank-N framework and
    **not** a Lambert-vs-Marshak principle; or

(b) the 1998 paper failed or was superseded, and Bogado Leite
    dropped the IC thread.

Either way, we learn that *in Bogado Leite's own career the 1998
paper is a terminus, not a springboard.*

### 3.2 Rodrigues 2011 thesis (OA, hdl.handle.net/10183/28965)

Retrieved the full PDF (8.6 MB, UFRGS `lume` repository) and
extracted via `pdftotext`. The thesis is primarily on direct grey
Dancoff calculations in cluster geometry via the PIJM algorithm
co-developed with Bogado Leite. Scan results:

- Bogado Leite is listed as co-advisor (`Co-orientador`).
- Reference list cites Bogado Leite 2000 / 2001 / 2003 / 2004a / 2004b.
- **Bogado Leite 1998 is NOT cited anywhere in the bibliography.**
- The "revised" interface-current relations are not referenced,
  not discussed, not transcribed. The thesis's theoretical core is
  scalar CP with Bickley $Ki_3$ integrals (Eq. 2.4 in the thesis)
  applied to cluster Dancoff factors — i.e. the classical Carlvik /
  Stamm'ler framework.

**Conclusion**: The 1998 revision has no living code descendant in
the Bogado-Leite-Vilhena research line. This is the strongest
ecosystem-level signal we can get without the PDF itself.

### 3.3 Sanchez — HAL archive scan

80-result pull on `authFullName_t:"Richard Sanchez"`; only one
paper in 2010–2016 with an OA document (2016 MOC tracking paper,
unrelated). **Sanchez 2014 NSE, his own 2014 paper, is not on HAL,
not on cea.hal.science, not on irsn.hal.science**. Standard Sanchez
has been actively depositing since the 2010s, so the 2014 NSE
paper's absence is a deliberate publication choice (likely at his
own discretion / ANS copyright policy at the time). This is a hard
block — there is no author-manuscript deposit to chase.

---

## 4. Sub-question verdicts

### Q1 — Does Sanchez 2014 give a rigorous IC-BC degeneracy theorem that makes the Lambert-Marshak asymmetry an allowed gauge choice?

**Partially addressed.** The Sanchez 2014 abstract establishes a
degeneracy of first-order $P_N$ equations *at the differential*
level and shows that IC/BC closures are a matter of principled
choice (via the second-order-parity equivalence). This is a
conceptual precedent for "multiple admissible IC/BC families exist
when the governing equations underdetermine them." But (i) F.4
operates on the **integral CP** formulation, not differential
$P_N$, and (ii) we do not have the main-text equations that would
let us map the solid-harmonic degeneracy argument onto the choice
of angular weight in $P_\mathrm{esc}$ vs $W$.

**Single verbatim excerpt useful to ORPHEUS**: the abstract's
closing sentence — "*our results reproduce those derived using
solid harmonic expansions by Davison and Rumyantsev in the 1950s*"
— tells us the underlying machinery is **solid harmonic
expansions**. Davison 1957 (*Neutron Transport Theory*, Oxford
Univ. Press) is in general circulation and would be the next
literature stop for anyone hunting the solid-harmonic derivation
Sanchez 2014 modernizes.

**Recommended action for Q1**: if the user can ILL the Sanchez 2014
PDF (ANS member library access is the cheapest route), we need the
main-text equations. The Davison 1957 monograph is likely a
simpler and open-access-available channel to the same underlying
solid-harmonic material.

---

### Q2 — Does Bogado Leite 1998 contain the exact Lambert-vs-Marshak derivation?

**Not in the literature checked.** Three independent lines of
evidence now weigh against this paper being a principled source:

1. **Zero abstract exposure** across CrossRef / OpenAlex /
   Semantic Scholar / Unpaywall (consistent with pre-abstract-
   deposition vintage, but unusual for ANE 1998).
2. **One citation total** on the paper (Semantic Scholar, still
   `citation_count: 1`, `influential: 0`), over 28 years.
3. **The author's own closest student (Rodrigues 2011 MSc) does
   not cite the paper**, and Bogado Leite's subsequent 2000–2004
   papers pivot to scalar CP + Dancoff with no IC follow-through.

The pre-Issue #122 memo
(`rank_n_ic_curvilinear_literature_leads.md`) already flagged
Bogado Leite 1998 as "the most dangerous false-negative." The
evidence marshalled in this pass confirms that while we cannot
*rule out* principled content, the ecosystem is silent. A PDF
recovery via ILL is the only path to certainty.

**Recommended action for Q2**: ILL request through the user's
institution, targeting the specific pages (ANE 25, issue 1-3,
pp. 129-139). CNEN/RJ (Rio de Janeiro) or IPEN/SP (São Paulo) may
also host an author-deposited PDF outside the public web index.

---

### Q3 — Does Corngold 2002 / 2004 addendum produce the 1/ρ prefactor of the diagonal $W_{oi,s}$ at $\sigma_t = 0$ from Peierls / Bickley-Naylor cylinder algebra?

**Not in the literature checked (at metadata level).** Both
Corngold papers are behind the Elsevier paywall with zero abstract
exposure. Title ("On the collision probability for the infinite
cylinder") and journal venue make it virtually certain that the
main content is the $Ki_n$ Peierls-integral analysis of the
cylinder CP. Whether Corngold's derivation produces a $1/\rho$
factor in the hollow-cell limit ($\sigma_t \to 0$, $r_\mathrm{inner}
= \rho R$) is an equation-level question that cannot be answered
without the PDF.

**Prior memory** (`phase4_cylinder_peierls.md` in this directory)
already contains the team's best public-literature synthesis of
cylinder Peierls algebra ($Ki_1$ vs $Ki_3$, $1/\pi$ prefactor,
chord branches). That memo is silent on the $1/\rho$ factor because
the public corpus there is silent on it. Corngold 2002 is the most
likely single source if the relation exists in print.

**Recommended action for Q3**: ILL the Corngold 2002 + 2004
addendum pair. The 2004 addendum is a short corrections note (cited
by zero), which suggests the body of the derivation is in the 2002
paper and the 2004 note fixes an algebraic error. Both are 6–15
pages each and are cheap PDFs.

---

## 5. Anything the literature supplies that the numerics-investigator can use today

This is the "even without PDFs, here is what we learned" section.
Three usable signals emerged from the metadata pass:

### 5.1 The degeneracy language is principled (Sanchez 2014)

Sanchez 2014 abstract is proof that the language of "multiple
admissible IC/BC families for an underdetermined differential
system" is an accepted mathematical framing in the reactor-physics
literature as of 2014. This means ORPHEUS's empirical finding
("Lambert P/G + Marshak W works while formally-consistent
alternatives don't") can be *framed* as a gauge choice without
ORPHEUS being the first to claim the framing. The Sanchez 2014
citation is a defensible literature hook for the ORPHEUS Sphinx
documentation of F.4.

**Action available today**: add a paragraph to ORPHEUS F.4 Sphinx
theory notes citing Sanchez 2014 abstract (not equations) as
precedent for "IC/BC closure choice in a degenerate governing
system may admit multiple principled forms." This is a publishable-
prose claim, not a symbolic derivation.

### 5.2 Davison 1957 monograph is probably the open-literature source of the solid-harmonic degeneracy argument

Sanchez's abstract explicitly *attributes* the underlying
derivation to Davison and Rumyantsev (1950s). Davison's *Neutron
Transport Theory* (Oxford 1957) is a well-known monograph and
should be in the user's Zotero (Zotero was flaky this session;
cannot confirm). If Davison's Chapter 6 or 7 contains the
solid-harmonic derivation of P_N interface conditions, that text is
the actual mathematical source and Sanchez 2014 is a synthesis.

**Action available today**: verify Davison 1957 Chapter 6 / 7 is in
Zotero. If yes, extract. If not, user should add — it is the
likely principled source Sanchez 2014 cites.

### 5.3 The 1998-vintage literature has convention-choice pathology documented empirically

Bogado Leite's 1999 TTSP abstract — "*differences in the Dancoff
factor approximations, in cylindrical and spherical geometries. …
significant discrepancies among the results*" — plus the title of
his 1998 paper ("**Revised** interface-current relations") is
independent evidence that *the particular convention choice in
cylindrical and spherical IC formulations was producing numerical
discrepancies across published formulations* by 1998. This is not
a principled derivation, but it is a *documented historical
precedent* for ORPHEUS's L2 observation.

**Action available today**: a future ORPHEUS Sphinx F.4 note may
cite Bogado Leite 1999 TTSP abstract as evidence that the
"convention-inconsistency causes significant discrepancies" phenomenon
was a live problem in the 1990s literature, even if we do not yet
have his 1998 equations.

---

## 6. Recommended next experiments (for the numerics-investigator, SymPy-ready)

The literature pass does not hand us a closed-form Lambert = Marshak
+ Δ identity. It does narrow what to test. Two candidate symbolic
experiments are ripe.

### 6.1 E-Q1: Solid-harmonic expansion of the outgoing angular flux on a sphere

**Motivation**: Sanchez 2014 abstract signposts that
Davison / Rumyantsev used *solid harmonic* expansions to derive
P_N interface conditions. Solid harmonics in 3D are
$Y_{\ell m}(\Omega)\, r^\ell$, with a natural split between even-ℓ
(symmetric under $\Omega \to -\Omega$) and odd-ℓ (antisymmetric).
The degeneracy Sanchez identifies is exactly that the first-order
P_N equations couple parity-mixed components that the integral form
treats independently.

**Test**: on the hollow-sphere outer surface, expand the outgoing
angular flux $\psi^+(\mu)$ (where $\mu = \Omega \cdot \hat n > 0$)
in two bases simultaneously:

- **Marshak basis**: shifted Legendre $\{P_n(2\mu - 1)\}$ on
  $[0,1]$, *with* µ-weight in the orthogonality measure. This is
  the P/G/W basis of formally-consistent rank-N closures.
- **Lambert basis**: also shifted Legendre $\{P_n(2\mu - 1)\}$ on
  $[0,1]$, but with *unweighted* orthogonality $\int_0^1 d\mu$. This
  is what F.4's P and G integrals use effectively (integrand
  $\sin\theta \cdot e^{-\tau}$, no explicit $\mu$ weight).

Compute the *change-of-basis* matrix $M_{nm} = \langle P_n^\mathrm{M},
P_m^\mathrm{L} \rangle$ and its inverse. The identity

$$
P_n^\mathrm{Marshak}(\mu) = \sum_m M_{nm}^{-1}\, P_m^\mathrm{Lambert}(\mu)
$$

is the **principled symbolic bridge** between the two conventions.
Whether F.4 can be read as "a basis rotation of formally-consistent
Marshak rank-N with an algebraic simplification that only closes
at rank-1" is something the numerics-investigator can check
analytically in SymPy in ~1 hour. **This was my single best yield
from the literature pass.**

### 6.2 E-Q2: Parity projection

**Motivation**: Sanchez 2014 uses "second-order equations with
opposite parity to N" as his equivalence construction. At rank-1
with $N = 1$, parity selects only even solid harmonics ($\ell = 0,
2, \dots$); at rank-2+ the odd-ℓ components start to contribute
independently to the interface conservation identity. This is a
*mechanism* for why rank-1 Lambert = Marshak (cancellation on
even-only modes) and rank-2+ Lambert ≠ Marshak (odd-ℓ modes leak
the gauge-freedom).

**Test**: symbolically evaluate $P_\mathrm{esc}$ and $W$ on the
hollow-sphere geometry in both conventions, project onto even-ℓ
and odd-ℓ solid harmonics separately, and check whether the
Lambert-Marshak **difference** lives entirely in the odd-ℓ
subspace at rank-1 (hence closes), and spills into the even-ℓ
subspace at rank-2 (hence opens the rank-2 plateau).

If this structural statement holds symbolically, it would explain
*both* L2 (why Lambert + Marshak works at N=1) and E2.4 (why the
naive Lambert-everywhere split-basis breaks). It would be a
negative-content result — no new closure, but a clean explanation
of the rank-1 specificity.

**Budget**: 2–3 hours SymPy, building on the existing
`diag_cin_aware_split_basis_keff.py` primitives.

### 6.3 Lower-priority backup: 1/ρ from Bickley $Ki_3$ at ρ → 0

**Motivation**: the L8 observation of a 1/ρ prefactor on the
diagonal $W_{oi,s}$ at $\sigma_t = 0$ is geometric, not algebraic,
and the Corngold 2002 derivation is the most likely public source.
Without Corngold's PDF, the user can reproduce the derivation from
scratch: in the hollow-sphere geometry with $r_\mathrm{inner}
= \rho R$, the inner-surface chord length averaged over isotropic
outgoing flux goes like $\bar \ell \propto 1/\rho$ as $\rho \to 0$
(small cavity ⇒ long mean chord on the inner surface per unit
inner area). The $Ki_3$ integrand $\sin^2\theta \exp(-\tau)$
limits to $\sin^2\theta$ at $\sigma_t = 0$, and the resulting
$\int_0^{\pi/2} \sin^2\theta\, d\theta \cdot \bar\ell$ gives the
$1/\rho$ factor. This is a standard Peierls calculation and does
not require Corngold. **Action**: hand the derivation to the
numerics-investigator as a scalar SymPy integral; a clean
symbolic proof is the deliverable.

---

## 7. Hand-off summary

**Addressed by public metadata**:
- Sanchez 2014 abstract (full) — gauge-theoretic precedent at the
  differential level, not the integral CP level. Q1 is
  partially addressed.
- Bogado Leite 1999 TTSP abstract — confirms
  convention-discrepancy phenomenology was documented by 1999.
  Q2 receives indirect corroboration but no principled derivation.

**Not addressed by public metadata** (PDFs blocked):
- Bogado Leite 1998 main text (Q2 direct).
- Corngold 2002 + 2004 main text (Q3).
- Wio 1984, Krishnani 1982, Krishnani 1985 main texts
  (transformation laws).

**Single highest-value follow-up**: ILL the Bogado Leite 1998 ANE
25(1-3) pp. 129-139 paper. Second-highest: Corngold 2002 ANE 29
+ 2004 addendum. Both via any AISU-member library (ANE is Elsevier;
any major research university has the subscription). If those
three PDFs arrive, this file should be updated with verbatim
equations and the Q1/Q2/Q3 verdicts upgraded.

**Single highest-value action NOT requiring any PDF**: solid-harmonic
change-of-basis between Lambert and Marshak conventions on the
outgoing-µ half-line [0, 1], as described in §6.1. SymPy, 1 hour,
yields a principled identity or a principled negative result.

---

## 8. Books and open references still to check (if the search continues)

These are the canonical monographs I expect to contain the
Davison / Rumyantsev solid-harmonic material Sanchez 2014
modernizes. I did not verify Zotero availability due to the
server's flaky state this session.

- **Davison, B. (1957)**. *Neutron Transport Theory*. Oxford
  University Press (also published as AEC document NAA-SR-3509).
  The standard English-language source for pre-1960 P_N interface
  conditions. Likely contains Chapter 7 / 8 solid-harmonic
  derivation. Out of copyright in some jurisdictions — OSTI may
  have the AEC report version.
- **Case & Zweifel (1967)**. *Linear Transport Theory*.
  Addison-Wesley. Chapter 2 covers analytical IC/BC constructions;
  Case-Zweifel completeness theorems are the direct ancestor of
  Sanchez 2014's degeneracy claim.
- **Williams, M. M. R. (1971)**. *Mathematical Methods in Particle
  Transport Theory*. Butterworths. The canonical bridge between
  differential P_N closures and integral CP.

**Open question for user**: are Davison 1957 and Case-Zweifel 1967
in Zotero? If not, both warrant addition. Davison's AEC report
version (NAA-SR-3509) may be OSTI-accessible as a 1957 AEC
document. This is a cheap Tier-2 follow-up (5 min) next session.
