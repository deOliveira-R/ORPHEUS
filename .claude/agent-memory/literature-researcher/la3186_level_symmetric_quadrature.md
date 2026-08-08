---
name: la3186-level-symmetric-quadrature
description: LA-3186 (Lathrop-Carlson 1964/65) + LA-4058 (Lathrop TWOTRAN 1968) full extraction for issue 337 — LS construction, axis-weight class ansatz, positivity frontier
metadata:
  type: project
---

# LA-3186 & LA-4058 — level-symmetric quadrature (issue #337, extracted 2026-08-08)

Full extraction: `scratch/issue_337_la3186_la4058_extraction.md` (page-cited,
render-verified). Reproduction script: `scratch/issue_337_cl_ansatz_check.py`.
Table dump: `scratch/la3186_tables_dump.txt`.

**Why:** #337 asked whether the S20/S22 positivity frontier is chaseable —
i.e. is the n > 12 point-weight class grouping systematic or ad-hoc table data.

**Identities.** LA-3186 = Lathrop & Carlson, "Discrete Ordinates Angular
Quadrature of the Neutron Transport Equation," written 1964-09-21, distr.
1965-02-12 (printed p = PDF p − 4). LA-4058 = Kaye D. Lathrop, "User's Guide
for the TWOTRAN (x,y) Program," 1968-12-01/1969-02-14 (printed = PDF − 3) —
a CONSUMER of LA-3186 (built-ins cap at ISN=16, quadrant Σw = 1/4, arrangement
cited to LA-3186 printed p. 10); it does NOT extend the LS tables.

**The construction (LA-3186)**: (a) level system Eq. (8) — even moments
i = 0..n/2 in (μ₁, level weights w) ONLY, μ₁ root-found on the top moment;
class structure never enters (8). (b) Point decomposition: p. 15 axis-weight
ansatz p_{orbit{i,j,k}} = a_i + a_j + a_k (i+j+k = n/2+2) — SYSTEMATIC, all n;
Fig. 3 (printed p. 10) instantiates classes ≤ n=16 with the relations
(n=14: p₃+2p₆ = p₄+p₅+p₇; n=16: p₄+p₅+p₈ = p₃+p₆+p₇) that are THEOREMS of the
ansatz. Table I tabulates n = 4,6,8,12,16,20 only; **n=20 p-column is BLANK**
(level weights only, no published point weights above n=16 anywhere).
Printed μ₁: S4 .3500212 / S6 .2666355 / S8 .2182179 / S12 .1672126 /
S16 .1389568 / S20 .1206033 — S6/S12/S16 have 1-ulp PRINT rounding slips
(correctly rounded: …54, …27, …69). Tables normalize Σ(level w) = 0.5 but
Σ_octant p = 1 (mixed convention, footnote-verified).

**Positivity (the frontier)**: "for n > 22 negative w_j" is about LEVEL
weights and verifies (min level w: +0.0016 at 22, −0.0169 at 24). Reproduction
[M]: ansatz POINT weights go NEGATIVE at n=18 (−0.0050), positive at n=20
(+0.00098), negative at n=22 (−0.0180); LP over the full decomposition kernel:
n=18 best min p +0.0078, n=20 +0.0048, **n=22 has NO nonnegative decomposition
(max-min −0.0027) — S22 intrinsically dead for the even-moment LS family; the
point-positive frontier is S20**. rank(orbit→level) = n/2−1 always; kernel
dim = #orbits−(n/2−1): 0 at n≤12, 1 at 14/16 (ansatz = printed member),
2/3/4 at 18/20/22.

**How to apply:** any future LS-quadrature question (constructions, μ₁ values,
positivity claims, LQ_n provenance) starts from the extraction file, not the
PDFs. LA-3186 Refs: method-of-moments general point-weight solve = Carlson
LA-3060 (1964); recursion source = Carlson-Lee LAMS-2573 (1961). Later
point-weight tables (S18+) live OUTSIDE these two reports (TWOTRAN-II
LA-4848-MS, DANTSYS-era constants) — unverified lead.
