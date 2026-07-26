---
name: dsa-primaries-full-extraction
description: "DSA issue #2 memo DONE (rev. 2, 7 sources): 5 primaries + Adams-Martin M4S (§7) + Warsa k-eigenvalue IRAM (§8) extracted to .claude/plans/dsa_literature_memo.md — Alcouffe print-sign errata, four-step + K_N table, gate formulas, NSE-147 twin-paper trap, absent papers"
metadata:
  type: project
---

**The DSA literature memo is COMPLETE, rev. 2** (2026-07-26):
`.claude/plans/dsa_literature_memo.md` — per-source extraction + synthesis for
issue #2 (consistent DSA, slab DD first). §§1–5 = the five primaries; **§7 =
Adams-Martin 1992 M4S; §8 = Warsa et al. 2004 k-eigenvalue IRAM** (both added
rev. 2). All equations spot-verified on rendered pages.

**Why:** Phase-3 DSA needs the exact discrete low-order operator, BC rows, and ρ
gates; the memo is the implementing session's sole literature surface.

**How to apply / durable traps (don't re-derive):**
- **Alcouffe 1977 PRINT SIGN SLIPS** (NSE 64:344): printed Eq. (17) carries spurious
  leading minuses and Eq. (23) prints "+R" where "−R" is required. Resolution pinned
  by continuous (2)-(3) [R = ∇·J̃ + ∇·D∇φ̃, SUBTRACTED], the convergence identities
  (19)/(20), the model-problem (24)-(26) [sign-consistent], and independently by
  Morel 1982 Eq. (20). Memo §1.5. Also: Alcouffe Eq. (16) coefficient is **⅓** —
  the sidecar OCR renders it ½ (the known ⅓→½ OCR hazard). Alcouffe gives NO
  boundary rows (deferred to ONETRAN-DA report, not local) and footnote-5 misprints
  Reed's year (1971, not 1973).
- **Larsen Part I** (NSE 82:47) = THE 3a reference: four-step = L₀+L₁ moments of
  balance AND closure (16a-d) → index promotion (18, lag φ₂ + closure residuals) →
  subtract, Schur-eliminate f₁ → tridiagonal (27) + updates (28); coefficients
  (23a-f) [DD: ρ=0]; BC rows (34)-(40) with discrete γ_N = Σ_{μ>0}μω (→¼, Σω=1);
  S2 ⟹ K₂=0 one-iteration exactness = the sharpest BC unit test; K_N Table I
  (WD .185→.248, LC ≤.315, LD ≤.300, LM ≤.282; K_N = 1−3ρ² closed form, K_∞ = ¼).
- **Normalization split**: Alcouffe/Larsen/McCoy/Morel/**Adams-Martin** Σw = 1
  (A-M proof: (25a) half-range Σw(n·Ω)² = 1/6, (25b) α ≈ ¼); A&L 2002 Σw = 2
  (their γ_N → ½). Map once before comparing any rows.
- **McCoy-Larsen Part II** (NSE 82:64): Table II = partial-consistency negative
  control (WD transport + DD low-order incl. BCs → divergence at σ_Th ≳ 10 for ANY
  a ≠ 0); Table V = P0-vs-P1 anisotropy ladder; the 4-region water/iron shielding
  problem fully spec'd in memo §3.1; fixup = "anathema" study Tables VI-IX.
- **Morel 1982** (NSE 82:34): P1 acceleration (27)/(28); curvilinear low-order uses
  cell-centered areas (37a) — NO unconditional-stability proof covers curvilinear or
  nonlinear-removal variants (p. 39-40); the ONLY explicit ALBEDO low-order rows in
  the corpus (Appendix, (1−α)/(1+α) factor); "consistency necessary but NOT
  sufficient" (p. 37); S_n diffusion theory D = 1/[3(σ_tr + {0,β/r,2β/r})].
- **A&L 2002**: reuse `.claude/plans/phase_i_survey_adams_larsen_2002.md` (full).
  Primary gate formula = (3.65) ω = ω_SI − c(1−ω_SI)/(1−c+⅓Λ²).
- **Adams-Martin 1992 M4S** (NSE 111:145, memo §7): the modification is STEP 2 of
  the four-step — ALSO lag the first-moment SURFACE-JUMP term C(r_k) (φ₀/φ₁
  interface discontinuities) ⟹ exact discrete Fick collapse (32) f₁ = −D∇f₀ ⟹
  step 4 always eliminable, any dim/grid; two-path equivalence (direct DFE of
  diffusion (18) ≡ four-step output (34)); slab-LD instance (35a/b/c) + BOTH-moment
  update (36); tridiagonal via App. C; spherical App. A ((A.6)+(A.9)) = the corpus's
  only curvilinear DFE-DSA operator. ρ: 0.2247c fine (0.0 S2) / peak 0.67c(S2)→
  0.50c(S16) at 1–2 mfp / **→0 thick** (unique). Paper says "not strictly
  consistent"/"inconsistent" — NEVER "partially consistent" (that label = A&L 2002
  §IV.D). DIFFERENT operator from Larsen §V fully-consistent LD (max 0.25 coarse);
  θ=1 LD variant of Larsen-Morel 1989 JCP 83 (the ORPHEUS slab-LD ref) behaves
  identically. Ref. 19 = Adams-Larsen 1989 Santa Fe = named LINEAR k-DSA primary.
- **⚠ THE NSE-147 TWIN-PAPER TRAP**: TWO "Warsa et al. 2004" NSE 147 papers.
  **147(1):26–42** = k-eigenvalue IRAM (5 authors, +McGhee+Lehoucq, DOI
  10.13182/NSE04-1) — LOCAL, memo §8. **147(2):218–248** = degraded-DSA robustness
  (3 authors) — STILL ABSENT (cited inside 147(1) as Ref [25] in submitted form,
  early title "DSA—Part I: Deficiencies in Multi-Dimensional Heterogeneous
  Problems"). Never let "Warsa 2004" resolve unchecked.
- **Warsa k-eigenvalue paper content** (memo §8): Arnoldi/IRAM (ARPACK) on
  **A = DH⁻¹F, H = L − MSD** (moment-vector eigenproblem); DSA = within-group WLA
  partially-consistent (CFEM diffusion by CG → DFEM update), in BOTH postures
  (PI: SI+DSA-acceleration+warm-start; IRAM: DSA-preconditioned FGMRES/BiCGStab +
  BFG relaxed-tolerance cascade, constant 0.1 — Arnoldi forfeits warm starts);
  IRAM outers FLAT in dominance ratio (14–26 @ δ=0.970–0.999 vs PI 468–8824,
  Table IV); MOX C5G7 28M DOF: 5× cheaper without upscatter, only ~75% with
  (nested upscatter = the declared drawback).
- **ABSENT from local library** (user to add if wanted): Warsa-Wareing-Morel NSE
  **147(2)**:218 (degraded DSA); Yavuz-Larsen 1988 Trans. ANS 56:305 (reflective-BC
  DSA).
- **Sidecar state**: #6/#7 (A-M, Warsa) have pypdf `.textlayer.md` files, NOT
  Mistral sidecars — the Mistral key went 401 mid-day 2026-07-26 (worked at 00:50).
  Failover that worked: pypdf text-layer extraction (same `## p. N` convention) for
  grep/navigation + `Read pages=` for verification. Regenerate true sidecars when
  the key is fixed (`ocr_literature.py --glob`, cache-idempotent).
