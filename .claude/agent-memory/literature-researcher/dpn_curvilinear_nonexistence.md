---
name: DP_N tables are slab-only — no cylinder/sphere analogue exists
description: Record the hard-searched negative result that DP_N k_eff reference tables (Stepanek-style) do not exist for curvilinear geometries; name the correct surrogates.
type: reference
---

**Fact**: DP_N / double-spherical-harmonics / Yvon method is a **slab-only**
construct — the double-Legendre basis splits at μ = 0 (inward/outward
half-ranges), which is meaningful only at a planar interface. No
published DP_N k_eff tables exist for solid cylinder or solid sphere
with white/reflective BC.

Checked exhaustively (2026-04-18 search): Stepanek 1981
(DOI 10.13182/NSE81-A19606, slab only), Stepanek 1986 TRISM follow-up
(DOI 10.1016/0149-1970(86)90023-5), Sanchez & McCormick 1982 review
(DOI 10.13182/NSE80-04-481 — DP_N appears only in slab context; §III.F
does NOT cover DP_N for 3-D geometries), Ganapol 2024 DP_N benchmark
(DOI 10.3390/foundations4030027, arXiv 2406.12942 — explicitly plane
geometry), Carlvik 1967 (CP only), Sood 2003 (F_N analytical).

**Correct surrogates** for DP_N verification in curvilinear geometry:

- **Sood, Forster, Parsons 2003** (Zotero F3D67AMQ,
  DOI 10.1016/s0149-1970(02)00098-7): 9 cylinder + 14 sphere one-group
  analytic k_eff to ≥5 decimals (F_N / singular eigenfunction), c =
  1.02–1.50. This is the community's canonical exact-truth set.
- **Ganapol 2024 Foundations** (DOI 10.3390/foundations4030027, OA):
  template methodology for tabulating DP_N convergence toward an F_N
  reference, in slab geometry. Port the methodology; port-able open
  reference implementation.

**How to apply**: when asked for "canonical DP_N tables" in a curvilinear
solver verification context, state the non-existence plainly and steer
toward (a) Sood 2003 exact values as truth and (b) locally-generated
DP_N truncation studies that demonstrate convergence toward those
values. Do not search further — the answer is settled.
