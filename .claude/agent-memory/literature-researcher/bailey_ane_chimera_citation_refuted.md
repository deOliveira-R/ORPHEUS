---
name: bailey-ane-chimera-citation-refuted
description: "Bailey/Adams/Yang/Zika (2009) Ann. Nucl. Energy 35, 1929-1936 is a CHIMERA — not a real record; every field traces to a different real paper. The equations it was cited for (Eq. 50 alpha recursion, Eq. 52 eta-ascending level order) are BMC 2010 NSE 165, which is LOCAL and already in refs.bib. Also carries the corrected cylindrical (mu,eta,xi) notation map."
metadata:
  type: project
---

# "Bailey … 2009 Ann. Nucl. Energy 35, 1929-1936" — REFUTED, it is a chimera

Full audit: `scratch/p4_bailey_citation_audit.md` (2026-08-27). This file is the
negative-lookup row + the corrected notation map.

## ⛔ The verdict — three real records fused into one fake one

| the fake entry's field | the REAL record it came from |
|---|---|
| title *"A piecewise linear discontinuous FE spatial discretization of the transport equation"* | **Bailey, Adams, Chang (2008)**, *"…of the Transport Equation **in 2D Cylindrical Geometry**"*, **LLNL-CONF-407632**, CONFERENCE, OSTI **952424**, no DOI, free at `osti.gov/servlets/purl/952424`. ⚠ authors are **Chang**, not Yang/Zika; the geometry clause was dropped. |
| authors Bailey, Adams, **Yang, Zika** | the separately-retracted **JCP 227:3738-3757 (2008)** *diffusion* paper `10.1016/j.jcp.2007.11.026` |
| **Ann. Nucl. Energy 35, 1929-1936** | **Zio & Zoia (2008)**, *"Bayesian inference of BWR model parameters by Markov chain Monte Carlo"*, `10.1016/j.anucene.2008.03.007` |
| year **2009** | none — all three are 2008; ANE vol. **35 = 2008** |

`[M]` denominator: CrossRef `journals/0306-4549/works?query.author=Bailey`
**total-results = 4** over ALL of ANE (Gallala 2017 v101, Ortega-Slaybaugh-
Brown-**Bailey**-Chang 2020 v138, Lujan 2024 v208, Smith 2025 v223) — **none a
finite-element paper, none at vol. 35**. T. S. Bailey has no PWLD paper in ANE.

## ✅ Where the equations REALLY live — all LOCAL, all already in `docs/refs.bib`

**Bailey, Morel & Chang (2010), NSE 165(2):149-169, `10.13182/NSE08-66`**
(`BaileyMorelChang2010`; local PDF; printed p. ≈ PDF p. + 147). §V = R-Z.

| eq | printed p. | what it is |
|---|---|---|
| **(11)/(12)** | 151 | SPHERE α + edge-cosine recursions, `−2μ_m w_m`, `Σw = 2`, seeds `α_{1/2}=0`, `μ_{1/2}=−1` |
| **(50)** | **156** | **R-Z α dome recursion** `α_{m+1/2,n} = α_{m−1/2,n} − μ_{m,n} w_{m,n}`, `α_{1/2,n}=α_{M+1/2,n}=0`. No factor 2 (`ΣΣw = 4π`, Eq. 51). |
| **(52)** | **157** | **R-Z per-ξ-level edge-cosine recursion** `μ_{m+1/2,n} = μ_{m−1/2,n} + w̄_{m,n}`, seed `−√(1−ξ_n²)`, terminus `+√(1−ξ_n²)`, weights normalised per level to `2√(1−ξ_n²)` ⟹ **THIS is the η-ascending per-level ordering**, not Eq. 50 |
| **(53)/(54)** | 157 | per-level WDD closure `ψ_{m,n} = τψ_{m+1/2,n} + (1−τ)ψ_{m−1/2,n}`; starting flux `ψ_{1/2,n}` obeys a Cartesian eq. (55) |
| **(74)** | **160** `[SCAN]` | **Morel-Montry τ for R-Z**: `τ_{m,n} = (μ_{m,n} − μ_{m−1/2,n})/(μ_{m+1/2,n} − μ_{m−1/2,n})`, "analogous to those defined by Eq. (42)" (the sphere's) |

⚠ **Eqs. (50) AND (52) are BOTH printed self-referentially — published typos**
(`α_{m+1/2,n} = α_{m+1/2,n} − …`). `[VERIFIED ON SCAN]` pp. 156/157: the
rendered pages carry them, so OCR is faithful and the JOURNAL is wrong. Correct
RHS subscript `m−1/2`, proven by the correctly-printed spherical twins (11)/(12).
**Three instances of the same typo in one paper** — see [[lessons-L010]].

## ⭐⭐ NOTATION CORRECTION — my own AGENT.md table is INVERTED for the cylinder

`[VERIFIED ON SCAN]` **Hébert 2009 Ch.3 (3.152)/(3.151)/(3.157), printed p. 91**
and **BMC 2010 Eq. (48), printed p. 156** use the **SAME** assignment:

> **μ = RADIAL · η = AZIMUTHAL · ξ = AXIAL (the level cosine)**

Hébert: `μ = √(1−ξ²) cos ω`, `η = √(1−ξ²) sin ω`; conservative form
`(μ/ρ)∂_ρ(ρφ) − (1/ρ)∂_ω(ηφ)`. BMC: `(μ/r)∂_r(rψ) − (1/r)∂_ω(ηψ) + ξ∂_zψ`,
`ξ = cos θ`, `n` = the ξ-level index.

⛔ **AGENT.md's table ("our `mu_x` = their η radial; our `mu_z` = their μ axial")
is BACKWARDS for both Hébert and Bailey.** So is ORPHEUS
`directional.py:292`'s comment attributing `(η,ξ,μ) = radial, azimuthal, axial`
to the "Bailey 2009 / Hébert convention" — a **second, independent defect**
found during this audit (ORPHEUS may name its own axes freely; the
*attribution* is false). `lessons.md` **L-001** is the correct one and always was.
⟹ Against Hébert/BMC: `ORPHEUS mu_x ≡ their μ`, `mu_y ≡ their η`,
`mu_z ≡ their ξ`.

## Unread lead (flagged to user, not acquired)

**LLNL-CONF-407632** (above) — Stone-Adams piecewise-linear DFEM on arbitrary
polygonal **RZ** meshes + thick-diffusion-limit asymptotics vs bi-linear DFEM.
⚠ It is a **SPATIAL** paper (L-005); the 9 citing sites all make **ANGULAR**
claims that BMC 2010 already owns and numbers. Not needed for the fix — but it
is the plausible origin of the stray title string, and it is relevant to the
curvilinear-LD thread (#158, [[issue-158-linear-discontinuous-cell-update]]).

See [[morel-montry-tau-angular-cell-edges]],
[[lathrop-2000-angular-scheme-taxonomy]], [[refs-bib-g2-corrections]],
[[sphere-sn-pole-closure-canonical]] (the entry-(B) retraction, L-003).
