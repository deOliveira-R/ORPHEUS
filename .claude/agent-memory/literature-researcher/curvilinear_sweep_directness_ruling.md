---
name: curvilinear-sweep-directness-ruling
description: TAXONOMY RULING — classical 1-D curvilinear (sphere/cyl) SN sweep is a DIRECT one-pass triangular solve; the starting-direction (pole) flux is NEVER seeded from the previous iterate in any standard published formulation. Per-source equation inventory (Hébert / BMC-2010 / Lathrop-2000 / Morel-1989) + the one documented lag locus (albedo/white outer BC, with the direct shooting alternative).
metadata:
  type: reference
---

# Curvilinear SN sweep directness — the taxonomy ruling (2026-07-01)

**Cite when**: any "is the SweepOperator a direct inverse?" architecture
question, seed-independence invariants, ORPHEUS "Carlson coupled-pole
seed" from previous-iterate ψ (ruled a MAL-FORMULATION, verdict (c)).

Companions: [[phase-d-carlson-coupled-pole]] (Hébert sphere §3.9.4 depth),
[[sphere-sn-pole-closure-canonical]] (citation-correction history),
[[morel-1989-si-vs-apply-equivalence]] (triangularity criterion).

## The ruling (all four sub-questions)

1. **Starting-direction equation is CLOSED** (redistribution-free
   streaming+collision), solvable by a spatial sweep within the same
   source iteration:
   - Sphere μ=−1: Hébert (3.432) continuous, (3.433)–(3.435) DD sweep;
     BMC-2010 Eq. (17) "slab equation" −ψ'_{1/2}+σ_tψ_{1/2}=½σ_sφ+½Q.
   - 1-D cylinder ω=π: Hébert (3.407) continuous ("angular
     redistribution term vanishes"), (3.408)–(3.410) DD sweep. Hébert
     attributes the zero-weight prescription (3.396) to **Alcouffe &
     O'Dell** [his ref 36].
   - R-Z per ξ-level: BMC Eq. (55) "Cartesian equation" with
     μ_{1/2,n}=−√(1−ξ_n²) (edge-cosine recursion Eq. 52).
2. **The full sweep is a direct one-pass triangular solve.** α seed is
   algebraic (α_{1/2}=0; α_{M+1/2}=0 emerges automatically — BMC
   Eq. 11). ψ seed comes from the starting-direction sweep. Each
   finite-weight ordinate reads only already-computed values: Hébert
   (3.436)/(3.438) numerators read φ_{n−1/2,i} (previous ordinate, SAME
   iteration); edge recursions (3.437)/(3.439)/BMC (16). Lathrop-2000
   p.247: "**sequential solution** of the equations"; the μ_{1/2}=−1
   coefficient (1−μ²_{m−1/2})=0 "results in the … **decoupling** of the
   equation set from the μ=−1 equation" (one-way, lower-triangular).
   Morel-1989 p.75: standard SN angular leakage matrix IS lower
   triangular → SI sweep = the system solve; only a FULL (Galerkin)
   matrix forces "off-diagonals on RHS and iterated on" with ρ→1 risk.
3. **Morel–Montry closure = same sequential recursion, different τ.**
   WDD ψ_m=τ_mψ_{m+1/2}+(1−τ_m)ψ_{m−1/2} (BMC 15); M-M weight
   **τ_m=(μ_m−μ_{m−1/2})/(μ_{m+1/2}−μ_{m−1/2})** (BMC 42) ⇔
   μ_m=τ_mμ_{m+1/2}+(1−τ_m)μ_{m−1/2} (43), from forcing M&M's β≡
   Σμ_m[α_{m+1/2}μ_{m+1/2}−α_{m−1/2}μ_{m−1/2}]=0 (41). NO simultaneous
   block, NO previous iterate. Table I: M-M sum = 0 to round-off.
   M&M primary = TTSP 13(5):615 (1984), **DOI 10.1080/00411458408211661**
   (CrossRef-verified; NOT in scratch/literature/).
4. **No published formulation lags the pole/starting flux.** Lathrop
   2000 (DOI 10.13182/NSE00-A2114, NSE 134:239–264) compares 5 schemes —
   diamond, M-M WDD, linear continuous (= "the original S_N", Carlson &
   Lathrop 1968 Ch.3 Gordon & Breach), linear discontinuous
   (Walters-Morel 1991), quadratic continuous — "All of the continuous
   approximations were initialized using the plane-geometry form of the
   transport equation that results when μ=−1" (p.251); none lags.
   The ONE standard lag locus is the **albedo/white OUTER BC** (Hébert
   3.415–3.416: J⁻ from previous iteration's J⁺), and even there Hébert
   gives the DIRECT alternative — the **shooting method** (3.417),
   2 sweep pairs + linearity. BC lag ≠ pole lag.

**Verdict (c)**: seeding the μ=−1 flux from the previous iterate's ψ is
an ORPHEUS mal-formulation, not a documented approximation. At the SI
fixed point it coincides (if the lagged ψ was itself properly swept),
but it breaks seed-independence/direct-inverse, changes the iteration
operator, and de-linearizes the sweep for Krylov wrapping. The standard
treatment has NOTHING to lag — the starting equation is closed by
construction ((1−μ²)→0), per-iteration cost O(ng·nx).

## α-recursion normalization hazard (extends L-001)

| Source | Recursion | Weights | Redistribution divisor |
|---|---|---|---|
| BMC sphere (11) | α_{m+1/2}=α_{m−1/2}−2μ_m w_m | Σw=2 | r(αψ−αψ)/w_m |
| BMC R-Z (50) | α_{m+1/2,n}=α_{m−1/2,n}−μ_{m,n}w_{m,n} | ΣΣw=4π | /(r·w_{m,n}) |
| Hébert sphere (3.424) | α_{n+1/2}=α_{n−1/2}−2𝒲_nμ_n | GL on [−1,1] | ΔS_i/(2𝒲_n) |
| Hébert cyl (3.399) | α_{p,q+1/2}=α_{p,q−1/2}**+**𝒲_{p,q}μ_{p,q} | per-level | **−**(1/𝒲_{p,q})[…] in (3.400) |
| Carlson/L&M convention (user's form) | α_{m+1/2}=α_{m−1/2}−μ_m w_m | Σw=2 | (ΔA/w_m)[…] |
| Larsen-Morel-2010 review sphere (1.23b), symbol **β** | β_{n+1/2}=−2Σ_{n'≤n}μ_{n'}w_{n'} ⇔ β_{n+1/2}=β_{n−1/2}−2μ_n w_n | Σw=2 | +(r/w_n)[βψ−βψ] in r²-form (1.23a); ≡ BMC sphere α exactly; ⚠ their α = SPATIAL WD weight (1.30) |

Factor-of-2 and sign move between the recursion and the balance-equation
divisor; α_{1/2}=0 and lower-triangularity are invariant across ALL.

## Local PDFs used (scratch/literature/)

- `Hebert(2009)Chapter3.pdf` pdf-pp.72–78 = book pp.138–144 (§3.9.3–3.9.4).
- `Bailey-Morel-Chang(2010)…pdf` pdf-pp.4–10 = journal pp.151–157
  (NSE 165:149–169, DOI 10.13182/NSE08-66 — CrossRef-verified).
- `Lathropp(2000)…pdf` (NSE 134:239–264) — abstract, pp.245–251 (§III–IV),
  p.264 refs.
- `Morel(1989)…pdf` p.75 (lower-triangular passage re-verified this session).

Lewis & Miller and Morel-Montry-1984 remain ABSENT from the folder;
prior ruling "L&M §4.5 not needed (Hébert identical)" re-affirmed; M&M
primary quotes still unavailable locally (BMC 41–42 + Lathrop p.250
restate its prescription — Morel is a BMC co-author).
