---
name: Sood-Forster-Parsons 2003 vs LA-13511 1999 — META-extraction comparison
description: Verdict on whether Sood 2003 journal version contains F_N / multi-region / anisotropic METHOD specs we currently believe must come from cited references. Answer: NO. Both papers are equation-DEFINITION + RESULTS-TABULATION; method derivations live in cited literature in BOTH versions. Use 1999 for tables, 2003 for an updated bibliography pointer + general multigroup k_inf appendix.
type: project
---

# Sood 2003 (Prog. Nucl. Energy 42:1, pp.55-106) vs LA-13511 (1999) — full PDF read

## Bottom-line verdict

**Sood 2003 is a condensed, lightly-updated journal port of LA-13511.**
It does NOT contain the F_N / Westfall-Metcalf / Wiener-Hopf / X-function / multi-region / anisotropic-dispersion METHOD machinery we suspected might be hiding in the journal version. Both papers stop at the same point: defining the transport equation per (group-count × geometry × scattering-order) family, then citing literature for the analytic technique. The implementer still needs the cited primary sources (KLL 1974, Grandjean-Siewert 1979, Siewert-Thomas 1986, Westfall-Metcalf 1972/1973, Sanchez-Ganapol 1983, Burkart-Ishiguro-Siewert 1976, Atalay 1996/1997, Bosler 1972, Ishiguro-Garcia 1978).

The 2003 journal text and the 1999 report are structurally **isomorphic**: same Section V/3 transport-equation overview, same Appendix A k_inf derivations, same 75-problem catalogue, same reference values. No new equations, no new method exposition, no Wiener-Hopf factorization, no F_N matrix construction, no multi-region X-function, no anisotropic dispersion relation Λ(z).

## Structural delta (2003 − 1999)

| Aspect | 1999 (LA-13511) | 2003 (Prog. Nucl. Energy) | Delta |
|---|---|---|---|
| Pages | 57 | 52 | -5 (compressed format) |
| Numbered equations (body) | ~32 (slab+two-group+matrix) | ~12 (1-12) | -20 (renumbering) |
| Appendix A equations | A.1-A.32 | A.1-A.59 | +27 (new General Multigroup section V) |
| References | 45 | 54 | +9 (V&V conference papers; Zoldi 2002 thesis; Lewis-Miller 1993) |
| Sections (top-level) | 12 (I-XII) | 8 (1-8 + Refs + App A) | renumbered |
| Test problems | 75 | 75 | 0 |

**The single methodologically meaningful 2003 addition** is Appendix A § V "General Multigroup Infinite Medium k_inf" (Eqs. A.55-A.59), credited to Larsen 1998 private communication. This gives the matrix form `S φ = (1/k) χ ν Σf φ` and the scalar k = (matrix inversion) form. **No spatial method content** — purely the infinite-medium k_inf algebra extended to G groups.

Everything else is renumbering or stylistic prose.

## Method specifications the implementer needs — where they actually live

| Need | 1999 reference | 2003 reference | In Sood paper? |
|---|---|---|---|
| Slab F_N derivation | [7] Grandjean-Siewert 1979 + [26] KLL 1974 | [14] + [35] | **NO — citation only** |
| Sphere F_N derivation | [26] KLL 1974 + [8] Siewert-Thomas 1986 | [35] + [15] | **NO — citation only** |
| Cylinder Westfall-Metcalf | [27]-[29] Westfall + Metcalf 1972/1973/1983 | [36]-[38] | **NO — citation only** |
| Multi-region/reflected slab | [33]-[37] Stewart-Metcalf, Forster, Ishiguro-Garcia | [42]-[46] | **NO — citation only** |
| Anisotropic-scattering slab | [13]-[15] Atalay, Dahl-Sjostrand; [22]-[24] Burkart, Boffi, Bosler | [20]-[22], [29]-[31] | **NO — citation only** |
| Two-group method | [37] Ishiguro-Garcia 1978; [33]-[36] Stewart, Forster | [46] + [42]-[45] | **NO — citation only** |
| Cylinder anisotropic | [31] Sanchez-Ganapol 1983 | [40] | **NO — citation only** |

**Verdict on the user's hypothesis**: NO, the method is not "actually there" in Sood 2003. Both papers are catalog + transport-equation-definition + results, NOT method exposition. Neshat-Maiorino is not cited in either version.

## Errata corrected?

The 1999 Eq. 28 numerator (`x_1 (ν_2 Σ_{2f} Σ_{21,s} + Σ_2^rem ν_1 Σ_{1f}) + x_2 (ν_1 Σ_{1f} Σ_{12,s} + Σ_1^rem ν_2 Σ_{2f})` over `Σ_1^rem Σ_2^rem - Σ_{12,s} Σ_{21,s}`) appears verbatim as 2003 Eq. A.11 — same algebraic form, same subscripts, same denominator. **No correction**: if Phase A flagged a typo, it persisted into the journal version unchanged.

One vestigial typo *introduced* by the 2003 renumbering: line 2322 still reads "to obtain the flux ratio, equations 23 and 24 are added", but 2003 has no equations 23 or 24 (the appendix uses `A.6`, `A.7`). This is a renumbering bug — copy-edit failure during journal port, not a math error.

## Cases delta

Identical 75-problem catalogue. Same naming convention (`PUa-1-0-IN`, etc.), same critical radii, same flux ratios, same precision (5+ decimals). The 2003 paper preserves higher-precision literature values (often 6-7 decimals) per the abstract claim.

## Recommendation for the implementer

**Use 1999 (LA-13511) as the primary catalogue** — it has marginally larger fonts, cleaner equation typesetting, and one extra Appendix-A equation count from the spelled-out matrix derivation. **Use 2003 only for** (a) the General Multigroup k_inf derivation in Appendix A § V (if a >3-group infinite-medium check is ever wanted), and (b) the more recent reference list pointing at Lewis-Miller 1993 ([50]) and Zoldi 2002 V&V thesis ([8]).

**Do NOT pivot to 2003 as a method source.** The KLL 1974 wall (Branch 1 via Case singular eigenfunctions, per existing memory `kaper_lindeman_leaf_1974_fn_method.md`) remains. Dropping Neshat-Maiorino 1980 was correct: Sood does not cite it and does not provide its content.

The user's suspicion was sound to test, but the answer is conservative: Sood 2003 confirms what 1999 already established — it's a TEST SET, not a method paper. Method machinery must be sourced from primary literature in BOTH cases.
