---
name: diffusion-limit-consistency-per-scheme-booleans
description: SHARP per-scheme/closure True/False verdicts for ORPHEUS issue #236 diffusion_limit_consistent trait, each pinned to an exact equation read from the LOCAL source PDF. Resolves the separability memo's ambiguous "(DD without fixup)" prose. SPATIAL axis = LMM-1987 JCP 69 Table I + Larsen-Morel 1989 Part II §IV; ANGULAR axis = Bailey-Morel-Chang 2010 NSE 165 Eqs.(41)/(42). Joint predicate = spatial AND angular; Cartesian collapses the angular factor (β≡0, no α-term).
metadata:
  type: reference
---

# Diffusion-limit-consistency per-scheme booleans (issue #236)

**Cite when**: setting/justifying the `diffusion_limit_consistent: bool`
trait on an ORPHEUS SN spatial scheme (`DiamondDifference`,
`LinearDiscontinuous`, `Step`) or angular closure
(`MorelMontryAngularSweep`, `IdentityAngularClosure`); building the
`diffusion_limit_consistent(scheme, closure)` pairing-validity
predicate. Companion to (and SHARPENS) [[space-angle-discretization-separability]],
whose Q4 prose "(DD without fixup)" left DD ambiguous. All verdicts
below were read VERBATIM from the local source PDFs — not inferred.

## The governing definition (what "consistent" means)

A spatial scheme's thick-cell asymptotic limit must be **a stable and
consistent discretization of the diffusion equation for the
leading-order scalar flux** (LMM-1987 §I, p.286: "if application of
this procedure yields a legitimate discretized diffusion equation,
then good numerical accuracy can be expected"). BMC-2010 adds the
**first-order** refinement: the diffusion approximation is itself
correct to first order (their Eq. 8, O(ε²) error), so FULL consistency
also needs the FIRST-order scalar flux to satisfy the correct diffusion
discretization (BMC p.2, p.151). Leading-order (ε⁰) and first-order
(ε¹) are SEPARATE bars — DD clears one but not the other on the angular
axis; the spatial schemes split on leading-order alone.

## SPATIAL axis — Larsen-Morel-Miller 1987 + Larsen-Morel 1989 Part II

Source: `scratch/literature/Larsen-Morel-Miller(1987)...pdf` (JCP
**69(2):283-324**, DOI 10.1016/0021-9991(87)90170-7) and
`scratch/literature/Larsen-Morel(1989)...II.pdf` (JCP **83(2):212-236**,
DOI 10.1016/0021-9991(89)90229-5). LMM-1987 **Table I** (journal p.287,
PDF page 4) is the canonical per-scheme verdict table. Each scheme has
**four** limits (cell-edge × cell-average) × (thick × intermediate).

LMM-1987 Table I reconstructed (verified from PDF page 4 text +
section bodies, OCR-corrected):

| Scheme        | Thick-Edge | Thick-Avg | Interm-Edge | Interm-Avg |
|---------------|-----------|-----------|-------------|------------|
| Diamond (DD)  | maybe(q=0)| **yes**   | **yes**     | **yes**    |
| Step          | yes(a)    | maybe(b)  | **no**      | maybe(b)   |
| Lund-Wilson   | no        | maybe(c)  | no          | no         |
| Castor        | no        | maybe(c)  | no          | no         |
| New (=LD-family)| **yes** | **yes**   | **yes**     | **yes**    |

Qualifiers (LMM-1987 p.287): (a) yes only if boundary sources
isotropic; (b) yes only if σ_a=Q=0 and (σ_t·h) constant; (c) yes only
if (σ_t·h) constant.

### (a) DiamondDifference (DD) → **True** (leading-order)

VERDICT: **True**. Bare DD's thick-limit cell-AVERAGE flux IS a valid
leading-order diffusion discretization. Pinned to **LMM-1987 Eq. (4.24)**
(journal p.298): `ψ_{m,j} = ½(Φ_{j+1/2}+Φ_{j-1/2}) + O(ε)` where
Φ_{j+1/2} satisfies the edge-differenced diffusion equation **Eq. (4.22)**
— "the diamond-differenced, cell-average fluxes have the thick diffusion
limit" (p.299, top). And in the INTERMEDIATE regime DD gets BOTH
edge+average: **Eq. (4.33)** is the correct diffusion equation, **Eq.
(4.34)** gives `ψ_{m,j}=φ⁽⁰⁾(x_j)+O(ε)`, `ψ_{m,j+1/2}=φ⁽⁰⁾(x_{j+1/2})+O(ε)`
→ "the diamond-differenced cell-edge and cell-average fluxes have the
intermediate diffusion limit" (p.300). So DD is in the GOOD category on
the **leading-order / cell-average** standard — Table I Diamond row is
"yes" for thick-Average, intermediate-Edge, intermediate-Average.

THE NUANCE that the "(DD without fixup)" prose conflated: the only DD
"maybe" is the THICK-regime CELL-EDGE flux, and it fails ONLY when
q_m≠0 (anisotropic incident boundary flux), producing cell-to-cell
oscillations (Eq. 4.23, p.298-299; "they do oscillate substantially if
these fluxes are highly anisotropic", p.299). This is the well-known DD
edge-oscillation, NOT a failure of the leading-order scalar-flux limit.
Under the issue-#236 definition ("thick limit is a consistent diffusion
discretization for the LEADING-ORDER SCALAR flux"), DD QUALIFIES → True.

FIRST-ORDER caveat (do NOT set a separate trait False on this account
for SPATIAL DD): LMM analyze the SPATIAL leading-order limit only. The
DD FIRST-order deficiency that BMC-2010 surfaces is an ANGULAR-axis
result (the curvilinear flux dip), NOT a spatial-DD verdict. On the
SPATIAL axis, DD's leading-order limit is valid → the spatial trait is
True. (If ORPHEUS ever wants a stricter spatial "first-order-consistent"
sub-trait, DD/Step would be False and only LD True there — but that is
NOT what #236's definition asks. Larsen-Morel 1989 Part II §III shows
DD's thick limit also has accurate diffusion BOUNDARY conditions, p.222;
DD is the reference scheme used for the "exact" solutions, p.309.)

### (b) Step → **False**

VERDICT: **False**. Step has **no intermediate diffusion limit** in
general. Pinned to **LMM-1987 Eq. (5.20)** and the conclusion sentence
immediately after (journal p.303): "These equations imply that the step
cell-average and cell-edge fluxes do NOT possess the intermediate
diffusion limit unless Q_j=σ_{a,j}=0, σ_{T,j}=constant, and h_j=constant."
The thick-regime cell-average is itself only "maybe" (Table I qualifier
b: yes only if σ_a=Q=0 and σ_t·h constant). For a GENERAL diffusive
problem (σ_a≠0, spatially varying cross sections) Step fails both
intermediate limits and the thick-average limit → **False**. This
matches the ORPHEUS docstring claim ("Step has no valid intermediate
diffusion limit while LD has all four"). Note Step's COMPENSATING
virtue per LMM §V: it is "extremely simple and always produc[es]
positive solutions" (p.301) — the positivity-vs-consistency trade.

### (c) Full LinearDiscontinuous (slope source carried) → **True**

VERDICT: **True** — LD possesses **all four diffusion limits**. Pinned
to **Larsen-Morel 1989 Part II §IV**, which analyzes LD BY NAME (p.215:
"the diamond difference (DD) and linear discontinuous (LD) spatial
differencing schemes"). The LD scheme is defined by **Eqs. (4.1)-(4.3)**
(two moments ψ̄_{m,j}, slope ψ̂_{m,j}; closure Eq. 4.3c `ψ̄∓ψ̂=ψ_{out}`;
θ free, θ=½ the natural finite-element value). The headline result is
**Eqs. (4.16)-(4.19)** (p.226-227): "**Equation (4.16a) is a stable and
consistently differenced version of the diffusion equation (1.3a), and
Eqs. (4.16b) and (4.16c) are very accurate representations of the
diffusion boundary conditions** (1.3b) and (1.3c)." LD therefore has
the thick diffusion limit AND the correct diffusion boundary
conditions — both #236 conditions met. (In LMM-1987 the precursor "New"
scheme, §VII Eqs. 7.1-7.4, is the one stated to "possess ALL FOUR
diffusion limits", p.308, Eq. 7.33/p.312 "this new differencing scheme
has all four diffusion limits"; Part II then names the production LD
family and proves the same.)

CRITICAL implementation gate (already in the ORPHEUS docstring,
`orpheus/sn/spatial/linear_discontinuous.py:12-28`): the "all four
limits" property belongs to **FULL LD with the slope SOURCE moment
Q̂=Σ_s·φ̂ carried** (Part II Eq. 4.2/4.3a-b feed Q_j AND Q̂_j into both
moment balances). The FLAT-source cut (Q̂=0) is O(h²) but NOT
diffusion-limit-consistent. Since the brief states the slope source is
**now threaded into the iterate**, the trait → **True** is correct for
the current ORPHEUS LD. If a flat-LD variant is ever reintroduced it
must carry `diffusion_limit_consistent=False`.

## ANGULAR axis — Bailey-Morel-Chang 2010 (LOCAL)

Source: `scratch/literature/Bailey-Morel-Chang(2010)...pdf` (NSE
**165(2):149-169**, DOI 10.13182/NSE08-66). Analyzes the SN equations
discretized **only in angle** (space kept continuous, p.150 step 1) →
proves the angular axis carries its OWN diffusion condition,
independent of the spatial axis. The relevant trait is the
**first-order β-condition** (their wording: "Full consistency in the
diffusion limit requires first-order asymptotic preservation", p.2).

Two-tier angular result (verified from PDF pages 5-7):
- **Leading-order (ε⁰)**: **Eq. (32)** is "the correct diffusion
  equation ... with ANY choice of weighting factors in the
  approximation of the angular flux" (p.153, after Eq. 32). So the
  ε⁰ limit is preserved by step, diamond, AND Morel-Montry alike.
- **First-order (ε¹)**: the second-order current **Eq. (40)** carries
  a contamination term; J⁽²⁾ satisfies Fick's law (→ correct
  first-order diffusion) iff the β-functional **Eq. (41)** vanishes:
  `Σ_{m=1}^M µ_m [α_{m+1/2} µ_{m+1/2} − α_{m-1/2} µ_{m-1/2}] = 0`
  ("This summation is the same expression that Morel and Montry called
  β and forced to zero", p.155). β≠0 for step/diamond → first-order
  error → the flux DIP at the origin.

### (a) MorelMontryAngularSweep → `beta_first_order_consistent = True`

VERDICT: **True**. Morel-Montry weights make β=0 EXACTLY. Pinned to
**BMC Eq. (42)**: `τ_m = (µ_m − µ_{m-1/2})/(µ_{m+1/2} − µ_{m-1/2})`,
"Morel and Montry have shown that Eq. (41) [β=0] is satisfied assuming
the following weighted diamond weighting factors" (p.155). Eq. (43)
`µ_m = τ_m µ_{m+1/2} + (1−τ_m) µ_{m-1/2}` is the equivalent
cell-center↔cell-edge cosine relation. Confirmed numerically by BMC
**Table I** (p.156): the M-M sum is "zero to round-off" (e.g. 8.33e-17
at S8) while step (~0.2) and diamond (slowly →0) do not vanish. So
M-M is the unique weighted-diamond closure that is first-order
diffusion-limit-consistent.

### (b) IdentityAngularClosure (Cartesian, no redistribution) → **True (vacuously)**

VERDICT: **True** — the angular β-condition is VACUOUS / trivially
satisfied. The user's reasoning ("no angular redistribution ⇒ β≡0") is
SOUND, and here is the exact mechanism: the β-functional Eq. (41) is
built ENTIRELY from the α coefficients (α_{m±1/2}), and the α's are the
**curvilinear angular-redistribution recursion** — they originate ONLY
from the `(1/r)∂_v(ηψ)`-type angular-derivative term in curved geometry
(BMC R-Z **Eqs. (49)-(50)**, p.156: the α_{m+1/2,n} appear from the
`(α_{m+1/2}ψ_{m+1/2} − α_{m-1/2}ψ_{m-1/2})/(r w_m)` redistribution term,
with the recursion `α_{m+1/2,n}=α_{m-1/2,n}−µ_{m,n}w_{m,n}`,
`α_{1/2,n}=α_{M+1/2,n}=0`). In **slab/Cartesian geometry there is NO
such term** — the streaming operator is plain `µ ∂_x`, so ALL α≡0 and
the Eq. (41) sum is identically zero term-by-term. There is no angular
edge flux to close, no weighted-diamond-in-angle, nothing to get wrong.
Leading-order ε⁰ is preserved "with any choice of weighting factors"
(Eq. 32) AND first-order ε¹ is automatic because β≡0. So
`IdentityAngularClosure` → **True** is correct: the angular
diffusion-limit condition is trivially met when there is no angular
redistribution. (This is the formal content of the
[[space-angle-discretization-separability]] Cartesian-collapse note.)

## JOINT predicate

```
diffusion_limit_consistent(scheme, closure)
    = spatial_condition(scheme) AND angular_condition(closure)
```

The diffusion limit is a property of the (spatial, angular) **PAIR**,
not either factor alone — this is the LMM(spatial)/BMC(angular)
two-paper split made into one predicate. The conditions **factorize**
(each provable on its own axis: LMM kills space with continuous angle;
BMC kills angle with continuous space) but are **jointly required** — a
good M-M angular closure paired with bare Step STILL breaks the limit
(spatial factor False), and full LD paired with step-in-angle on a
curvilinear mesh STILL dips (angular factor False). Independence of
SELECTION ≠ independence of CONSEQUENCE.

**Cartesian collapse**: `IdentityAngularClosure` ⇒
`angular_condition = True` (vacuous, β≡0 because all α≡0, no
redistribution term in `µ∂_x`). So in Cartesian the predicate reduces
to the spatial condition alone:
`diffusion_limit_consistent = spatial_condition(scheme)`. This is the
formal justification for collapsing the angular factor in Cartesian.

## The trait table (issue #236)

| Item | axis | trait | value | citation |
|------|------|-------|-------|----------|
| `DiamondDifference` | spatial | `diffusion_limit_consistent` | **True** | LMM-1987 Eq.(4.24)/(4.22) thick-avg + Eq.(4.33)/(4.34) intermediate; Table I "Diamond" row |
| `Step` | spatial | `diffusion_limit_consistent` | **False** | LMM-1987 Eq.(5.20) + p.303 "do NOT possess the intermediate diffusion limit unless Q=σ_a=0…"; Table I "Step" no/maybe |
| `LinearDiscontinuous` (full, slope source) | spatial | `diffusion_limit_consistent` | **True** | Larsen-Morel 1989 Part II §IV Eqs.(4.1)-(4.3), result Eqs.(4.16)-(4.19) p.226-227 "stable and consistently differenced version of the diffusion equation"; LMM-1987 §VII "all four diffusion limits" |
| `MorelMontryAngularSweep` | angular | `beta_first_order_consistent` | **True** | BMC-2010 Eq.(42) τ_m sets β-functional Eq.(41)=0; Table I sum ~0 to round-off |
| `IdentityAngularClosure` | angular | `beta_first_order_consistent` | **True** (vacuous) | BMC-2010 Eq.(41) built from α's; α≡0 in Cartesian `µ∂_x` (no redistribution; cf. R-Z Eqs.49-50 where α's arise) ⇒ β≡0 |

⚠ DO NOT conflate the SPATIAL DD verdict (True, leading-order, LMM)
with the ANGULAR DD-in-angle verdict (β≠0, first-order-FAILING, BMC).
They are different axes; DD is True on one and the loser on the other.
The flux dip is an ANGULAR artifact, not a spatial-DD failure.

## Caveats / OCR provenance

- LMM-1987 Table I OCR is garbled in the raw extract (page 4); the
  table above is reconstructed by cross-reading each section's
  conclusion (§IV DD p.298-300, §V Step p.303, §VI LW/Castor p.308,
  §VII New p.312). The qualifier footnotes (a/b/c) are transcribed
  verbatim from p.287.
- All equation numbers above were read from the LOCAL PDFs this
  session; no Tier-2 lookup was needed (all three papers local).
- Zotero NOT consulted this session (verdicts come straight from the
  primary PDFs, the authoritative record).

## Cross-refs

- [[space-angle-discretization-separability]] — the parent memo; this
  one SHARPENS its Q4 "(DD without fixup)" prose into booleans.
- [[spherical_sn_central_cell_spatial_order]] — BMC is ANGULAR not
  spatial; the spatial closures (Hébert 3.9.4 / Stacey 9.9).
- [[issue_158_linear_discontinuous_cell_update]] — slab LD = Larsen-Morel
  1989 Eqs.(4.1)-(4.3); the slope-source threading gate.
- [[phase_d_carlson_coupled_pole]] — the curvilinear α-cascade /
  starting-direction sweep that GENERATES the α's in BMC Eq.(41).
