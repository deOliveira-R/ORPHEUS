---
name: face-transmission-pade
description: SN cell face transmission as Padé stability function — LR 1974 Thm 2 (DG(k) = [k/k+1], strongly A-stable), MWS 1996 Eq.74 (LD = (6-2τ)/(6+4τ+τ²)), DD [1/1] sources, multi-D DD −1 eigenvalue unpublished
metadata:
  type: reference
---

Extracted 2026-09-23 for the face-transmission AOR; memo was in session scratchpad (lit_face_transmission.md).

- **Lesaint-Raviart 1974** — FREE Rennes seminar version, numdam `PSMIR_1974___S4_A8_0`, now LOCAL `scratch/literature/Lesaint-Raviart(1974)...pdf` (no OCR sidecar yet: Mistral 429; pdftotext works). §2 ODE DG: Lemma 1 (= implicit RK), Lemma 2 (2.24) R=P/Q deg k/k+1, def (2.27) "strongly A-stable" (= L-stable), **Theorem 2 printed p.12: order 2k+1, R = "subdiagonal (k+1,k) Padé"** — LR order is (den,num), MWS is (num,den). Gauss-Radau nodes with ξ₁=0; LR never say "Radau IIA". Book version (de Boor, pp.89-123) NOT seen.
- **MWS 1996 JCP 128:445, PDF p.8 = printed 452, Eqs. (72)-(76)** [rendered]: LD "Padé (1,2)", negative for τ>3; lumped LD "Padé (0,2)".
- DD [1/1]: Stacey (9.210) printed 359 + negativity Δ>2|μ|/Σ; Hébert (3.489) printed 153 (Padé via Suslov); LMM 1987 JCP 69 (4.11)-(4.12) printed 297: thick-limit edge anisotropy alternates (−1)^j.
- ⛔ No local/primary source names DD "trapezoidal/Crank-Nicolson" or step "backward Euler"; Hairer-Wanner not local.
- ⛔ Multi-D DD octant transmission T = 2·1wᵀ − I (eig −1 ×(d−1) exactly, one eig = [1/1] at τ_eff=1/Σ|Ω_i|/(σh_i)) — NOT found published; ORPHEUS's own.
- Response matrix: Sanchez-McCormick 1982 §III.F.1 printed p.523: R = P_SS + P_S·V[Σ(1−PVΣ)⁻¹P·S], first term uncollided.
