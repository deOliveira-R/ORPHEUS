---
name: Kaper-Lindeman-Leaf 1974 acquisition attempt + F_N method literature map
description: KLL NSE 54, 94 (1974) is paywalled with no preprint and no OA mirror. Sood [26] cites it for VALUES only; the F_N method derivation lives in Garcia-Siewert NSE 69 (Part I+II, 1979) — also paywalled. Maps the access cliff and the secondary-source landscape.
type: reference
---

## Bibliographic record (verified)

- **Authors**: H. G. Kaper, A. J. Lindeman, G. K. Leaf (all Argonne National Laboratory)
- **Title**: "Benchmark Values for the Slab and Sphere Criticality Problem in
  One-Group Neutron Transport Theory"
- **Journal**: Nuclear Science and Engineering, vol. 54, no. 1, pp. 94-99,
  May 1974
- **DOI**: `10.13182/nse54-94` (NOT `NSE74-A23308` — that guess was wrong)
- **OSTI ID**: 4303686 (bibliographic only — NO fulltext, NO preprint)
- **Length**: 6 pages (94-99) — short Technical Note format
- **Cited by**: 39 (OpenAlex), 35 (S2). Influential: 0.

## Acquisition attempt — outcome: FAILED

| Source                   | Result                                                 |
|--------------------------|--------------------------------------------------------|
| Tandfonline (NSE host)   | Cloudflare 403 — paywalled                             |
| OpenAlex `oa_url`        | None — `is_oa: False`, sole location is the publisher  |
| Semantic Scholar PDF     | None                                                   |
| OSTI fulltext            | 404 (`/servlets/purl/4303686`); record bibliographic-only |
| OSTI ANL preprint search | NO ANL-XXXX preprint exists (Kaper bibliography 1972-85 surveyed) |
| INIS                     | API down at lookup time (TLS reset)                    |
| HAL / Zenodo / arXiv     | Out of scope (US, 1974)                                |

OpenAlex `locations[]` lists exactly one location: the NSE landing page,
not OA. There is no institutional repository deposit, no Argonne Reports
mirror, no author homepage copy.

**Conclusion**: KLL 1974 itself is unobtainable through the open
literature pipeline.

## Critical re-read of LA-13511 (Sood-Forster-Parsons 1999)

The action item assumed KLL contains the F_N method derivation.
**It does not.** KLL is the source of the *numerical values* — not
of the method. Sood's reference apparatus separates the two:

- **Sood Ref [26]** (Kaper-Lindeman-Leaf 1974) — the VALUES Sood
  tabulates for `Ua-1-0-SL` and `Ua-1-0-SP`.
- **Sood Ref [7]** (Grandjean & Siewert, "The F_N Method in Neutron
  Transport Theory. Part II: Applications and Numerical Results", NSE
  69, 161 (1979) — DOI `10.13182/nse79-a20608`) — the METHOD.
- **Sood Ref [6]** (Case & Zweifel 1967, "Linear Transport Theory")
  — the singular-eigenfunction substrate.

This shifts the acquisition target. The F_N derivation we need is in
the **Garcia-Siewert F_N pair** (Part I = `10.13182/nse79-1`, Part II =
`10.13182/nse79-a20608`). Both are also paywalled with no OA mirror
indexed by OpenAlex / S2 / OSTI.

## Secondary-source landscape — what reproduces F_N

Surveyed local PDFs and citing literature:

1. **Sood-Forster-Parsons LA-13511 (local)**: tabular values + brief
   sentence ("The analytic methods used include Case's singular
   eigenfunction, F_N and B_N methods, and Green's functions"). NO
   F_N derivation. No use beyond cross-checking final tabulated
   numbers.

2. **Valougeorgis (1985) "The F_N Method in Kinetic Theory" (PhD
   thesis, Virginia Tech, 297 pages)** — DOWNLOADED to
   `/workspaces/ORPHEUS/scratch/literature/Valougeorgis(1985)The Fn
   method in kinetic theory PhD thesis.pdf` (4.2 MB).
   - Chapter 3 "Formulation of the F_N Method" gives the
     **structural skeleton**: full-range orthogonality (Kriese et al.),
     singular integral equations for boundary distribution functions,
     finite expansion in basis functions, collocation. Generic to
     transport.
   - **Caveat**: this is BGK linearized gas dynamics, not one-group
     neutron transport. Basis functions, dispersion relation, and
     half-range factors differ. Equations are not directly
     transcribable. But the **method skeleton is identical** —
     useful as a roadmap to mimic, not as a code reference.
   - Cites Siewert et al. [7,8] for the original NSE 1979 form.
   - Plane-geometry development in §3.3-§3.4. Cylindrical in §6.

3. **Garcia spherical PDFs (local 2019/2021)** — use **P_N**, NOT
   F_N. They cite a *different* Kaper paper (Kaper-Shultis-Veninga
   1970 on slab albedo). NOT useful for KLL benchmark reproduction.

4. **Bell & Glasstone (1970)** — predates F_N method (1979). Has
   Case eigenfunction theory and slab/sphere critical problem
   solutions via singular eigenfunctions, but not the F_N
   approximation that yields the precise digits Sood reports.

5. **Williams (1971) "Mathematical Methods in Particle Transport
   Theory"** — predates F_N. Same caveat as B&G.

6. **Case-Zweifel (1967)** — Sood Ref [6], the singular-eigenfunction
   foundation. Slab benchmarks via Case's method exist there but at
   lower precision than KLL's F_N values.

## Assessment — can we proceed without KLL?

**Yes, for VALUES.** Sood/LA-13511 already publishes KLL's
five-digit-or-better numbers verbatim. The user has those values in the
brief: `r_c = 0.93772556 mfp` (`Ua-1-0-SL`) and `r_c = 2.4248249802 mfp`
(`Ua-1-0-SP`), plus flux ratios. Reading KLL directly would only confirm
provenance, not provide more digits.

**Acquisition not strictly required for Branch 1 SymPy IF the goal is to
match the published numbers**. Branch 1 of `algebra-of-record` is
"reproduce the published artefact bit-for-bit" — for that you need the
*method* (F_N), not the original tabulation paper.

**Acquisition IS required IF the goal is to re-derive F_N matrix
entries symbolically.** That requires Garcia-Siewert NSE 69 Part I/II,
both also paywalled.

## Recommendations for the user

1. **Drop Branch 1 = "reproduce KLL via F_N derivation".** That branch
   is blocked by the paywall on Garcia-Siewert NSE 69 + KLL NSE 54.

2. **Pivot Branch 1 = "reproduce KLL values via Case singular-eigenfunction
   integral equation"**. The eigenvalue equation for the bare-critical
   slab and sphere has a closed form via Case's method (Case-Zweifel
   1967, Bell-Glasstone Ch. 2.6, Williams Ch. 6). For one-group
   isotropic scattering with `c = ν Σ_f / Σ_t`, the critical equation
   reduces to a transcendental relation in `R/ν_0` involving the
   Case dispersion function `Λ(ν)` and `tanh(R/ν_0)` (slab) or
   `tanh(R/ν_0) - R/ν_0` (sphere). SymPy + mpmath can solve these to
   match KLL's published digits without ever touching F_N. This is
   the path of least resistance.

3. **Try institutional access for KLL/Garcia-Siewert if available**.
   The user's institution may have ANS NSE access via Tandfonline.
   That's an out-of-band channel; not something this agent can do.

4. **Valougeorgis 1985 was downloaded** as a method-skeleton reference.
   Useful when the user wants to *understand* F_N abstractly. Not
   useful for transcribing neutron-transport-specific equations.

## Specific numerical values requested

KLL was not acquired, so the digits cannot be confirmed from primary
source. The values in the user's brief match Sood Table for `Ua-1-0-SL`
and `PUb-1-0-SP` — those are KLL's tabulated values reprinted.
**Trust Sood's transcription** (it is itself a peer-reviewed report
funded by LANL with cross-checked digits, OSTI 762839 fulltext mirror).

## Action items if pursued further

- Try ANL Reports library directly (`https://anl.gov/library`) for
  any 1974 internal report by Kaper-Lindeman-Leaf that may be the
  preprint of NSE 54-94. OSTI did not return one.
- Try Argonne National Laboratory Technical Reports series (ANL-XXXX,
  ANL-7910 era) for the same year.
- Email the corresponding author or successors. Kaper has an ORCID
  (0000-0003-1154-2406) but is retired; G. K. Leaf may also be ANL
  emeritus.
- ILL request through the user's library — KLL is short (6 pp) and
  trivial to scan-and-deliver.
