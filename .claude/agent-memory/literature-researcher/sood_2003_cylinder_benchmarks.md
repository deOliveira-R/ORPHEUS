---
name: Sood 2003 cylinder critical-radius benchmarks (LA-13511)
description: Bare-cylinder one-group r_c values to 6 digits for Variant α α=0 cross-check; structurally independent of Bickley-Naylor (F_N method).
type: reference
---

Sood, Forster, Parsons — LA-13511 (1999) full report + *Prog. Nucl. Energy* 42(1) 55-106 (2003),
DOI 10.1016/S0149-1970(02)00098-7. Open-access at OSTI 762839 (paper) and 10601 (LA-13511).

Local cached PDFs (in transient /tmp; re-download from OSTI as needed):
- https://www.osti.gov/servlets/purl/762839 (paper, 601 kB)
- https://www.osti.gov/servlets/purl/10601 (LA-13511 long, 2.5 MB)

Cylinder benchmarks (one-group, isotropic scattering, vacuum BC, 1-D radial, infinite axial):

LA-13511 Table 13 — `Ua-1-O-CY` (U-235, c=1.30):
  ν=2.70, Σ_f=0.065280, Σ_c=0.013056, Σ_s=0.248064, Σ_t=0.32640 cm⁻¹.
  r_c = **1.72500292 mfp = 5.284935 cm**. Refs 27, 28 (Westfall 1983 / Westfall-Metcalf 1972, F_N method).

LA-13511 Table 7-8 — `PUb-1-O-CY` (Pu-239 (b), c=1.40):
  r_c = 1.396979 mfp = 4.279960 cm. Flux ratios Φ(0.5 r_c)=0.8093, Φ(r_c)=0.2926 (4 dp only).

Cylinder flux profile NOT tabulated for U-235 case (Table 14 has slab + sphere only).

Underlying method = F_N (Case singular eigenfunctions / Wiener-Hopf). NO Bickley-Naylor Ki_n.
Structurally independent of ORPHEUS Variant α's Ki_3-based T_00 path → eliminates ERR-032 reuse risk.

ORPHEUS use: drive Variant α at α=0, sweep R, find k_eff(R)=1, assert R/mfp = 1.72500292 ± 1e-5.
V&V pillar: 2 (semi-analytical). 5-digit budget per Sood's stated accuracy.

Gap (worth recording as V&V info): NO external published one-group homogeneous cylinder
k_eff at α∈(0,1). Sood α=0 and V_α1 algebraic α=1 bookend the parameter range; interior is
testable only via in-codebase cross-code (S_N or MOC, Pillar 4) if those solvers support
specular cylinder BC with tunable α.
