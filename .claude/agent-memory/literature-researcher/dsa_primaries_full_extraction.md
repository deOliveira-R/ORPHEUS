---
name: dsa-primaries-full-extraction
description: "DSA issue #2 memo DONE: all 5 primaries (Alcouffe/Larsen-I/McCoy-Larsen-II/A&L-2002/Morel) extracted to .claude/plans/dsa_literature_memo.md — Alcouffe print-sign errata, four-step + K_N table, gate formulas, absent papers"
metadata:
  type: project
---

**The DSA literature memo is COMPLETE** (2026-07-26):
`.claude/plans/dsa_literature_memo.md` — per-source extraction + synthesis for
issue #2 (consistent DSA, slab DD first). All equations spot-verified on rendered
pages. All 5 primaries LOCAL with sidecars.

**Why:** Phase-3 DSA needs the exact discrete low-order operator, BC rows, and ρ
gates; the memo is the implementing session's sole literature surface.

**How to apply / durable traps (don't re-derive):**
- **Alcouffe 1977 PRINT SIGN SLIPS** (NSE 64:344): printed Eq. (17) carries spurious
  leading minuses and Eq. (23) prints "+R" where "−R" is required. Resolution pinned
  by continuous (2)-(3) [R = ∇·J̃ + ∇·D∇φ̃, SUBTRACTED], the convergence identities
  (19)/(20), the model-problem (24)-(26) [sign-consistent], and independently by
  Morel 1982 Eq. (20). Memo §1.5. Also: Alcouffe Eq. (16) coefficient is **⅓** —
  the sidecar OCR renders it ½ (the known ⅓→½ OCR hazard, now confirmed in TWO
  papers). Alcouffe gives NO boundary rows (deferred to ONETRAN-DA report, not
  local) and footnote-5 misprints Reed's year (1971, not 1973).
- **Larsen Part I** (NSE 82:47) = THE 3a reference: four-step = L₀+L₁ moments of
  balance AND closure (16a-d) → index promotion (18, lag φ₂ + closure residuals) →
  subtract, Schur-eliminate f₁ → tridiagonal (27) + updates (28); coefficients
  (23a-f) [DD: ρ=0]; BC rows (34)-(40) with discrete γ_N = Σ_{μ>0}μω (→¼, Σω=1);
  S2 ⟹ K₂=0 one-iteration exactness = the sharpest BC unit test; K_N Table I
  (WD .185→.248, LC ≤.315, LD ≤.300, LM ≤.282; K_N = 1−3ρ² closed form, K_∞ = ¼).
- **Normalization split**: Alcouffe/Larsen/McCoy/Morel Σw = 1; A&L 2002 Σw = 2
  (their γ_N → ½). Map once before comparing any rows.
- **McCoy-Larsen Part II** (NSE 82:64): Table II = partial-consistency negative
  control (WD transport + DD low-order incl. BCs → divergence at σ_Th ≳ 10 for ANY
  a ≠ 0); Table V = P0-vs-P1 anisotropy ladder; the 4-region water/iron shielding
  problem fully spec'd in memo §3.1; fixup = "anathema" study Tables VI-IX.
- **Morel 1982** (NSE 82:34): P1 acceleration (27)/(28); curvilinear low-order uses
  cell-centered areas (37a) — NO unconditional-stability proof covers curvilinear or
  nonlinear-removal variants (p. 39-40); the ONLY explicit ALBEDO low-order rows in
  the corpus (Appendix, (1−α)/(1+α) factor); "consistency necessary but NOT
  sufficient" (p. 37) — the operative property is derived-by-moment-reduction;
  S_n diffusion theory D = 1/[3(σ_tr + {0,β/r,2β/r})].
- **A&L 2002**: reuse `.claude/plans/phase_i_survey_adams_larsen_2002.md` (full).
  Reflective-BC DSA theory is NOT in the review body — Yavuz-Larsen 1988 Trans. ANS
  56:305 [135] is bibliography-only. Primary gate formula = (3.65)
  ω = ω_SI − c(1−ω_SI)/(1−c+⅓Λ²), Λ = (2/Σ_th)tan(Σ_thλ/2).
- **ABSENT from local library** (user to add if wanted): Adams & Martin 1992 NSE
  111:145 (M4S); Warsa-Wareing-Morel 2004 (Krylov-DSA robustness); Yavuz-Larsen
  1988 (reflective-BC DSA).
