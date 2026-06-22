---
name: peierls-mg-fission-source-chi-indexing
description: Canonical literature ground that the fission emission spectrum χ in the multigroup integral (Peierls) transport source is SOURCE-indexed (evaluated at the fissioning point r_j, same r as νΣ_f and φ) — for ERR-063 / ORPHEUS #257.
metadata:
  type: project
---

ERR-063 / ORPHEUS #257 (Peierls MG fission χ sink/source swap) settled
SOURCE-indexed (convention B correct, A wrong) on literature ground.

**Canonical reference**: Hébert (2009) *Applied Reactor Physics* Ch.3,
the integral-transport + source-density chapter. Local PDF:
`scratch/literature/Hebert(2009)Chapter3.pdf`.

- Integral (Peierls) form: **Eq. (3.40)** φ(r,Ω,t)=∫₀^∞ ds e^{−τ(s)} Q(r−sΩ,Ω,t−s/v)
  — the kernel transports the emission density Q evaluated AT THE SOURCE
  POINT r−sΩ. Q is one well-defined local emission at the source point.
- Isotropic fission source: **Eq. (3.57)** Q_fiss(r,E)=Σ_j χ_j(E) ∫dE' νΣ_{f,j}(r,E') φ(r,E').
  EVERY spatial argument is the SAME r: χ_j (isotope-j spectrum), νΣ_{f,j}(r,·),
  φ(r,·). Fission source is purely LOCAL at r.
- Complete multiplicative source: **Eq. (3.58)** — scatter (first-index=source,
  Σ_{s,ℓ}(r,E→E')) + (1/4πk_eff)Σ_j χ_j(E)∫dE' νΣ_{f,j}(r,E') φ(r,E').
- χ definition: text above Eq. (3.56): χ_i(E)dE = probability an emitted neutron
  has energy E for fissionable nuclide i; normalization ∫χ_i(E)dE=1 (3.56).
  χ is a PROPERTY OF THE FISSIONING NUCLIDE → indexed to the fission/source point.
- Adjoint (3.61): χ_j(E') sits with φ*(r,E') — arguments interchanged under
  adjoint (rule 3 p.77), confirming forward χ shares the source/fission index.

**Corroboration** (structurally independent): differential-form fission
source S_g(x)=(χ_g/4π)Σ_{g'}νΣ_{f,g'}Φ_{g'} is single-point — all factors carry
ONE x. In Peierls form that whole local point-quantity IS Q(r_j); kernel then
streams it to r_i. So χ canNOT carry the sink index r_i.

**Verdict**: χ is SOURCE-indexed. (A) χ_{g}(r_i)·Σ νΣ_f(r_j)φ(r_j) is WRONG;
(B) χ_{g}(r_j)·Σ νΣ_f(r_j)φ(r_j) is correct. Differs only when χ varies in
space (multi-region, differing emission spectra); coincides for uniform χ.

**Zotero was DOWN** this session (port 23119 conn-refused, MCP tools not
surfaced — the documented flakiness pattern). Ground is from the local PDF
corpus + Tier-2 corroboration; no user annotations checked. Bell & Glasstone
1970 (OSTI 4074688 / UNT mirror) is the foundational χ-as-nuclide-property
source but the UNT PDF is CAPTCHA-gated.
