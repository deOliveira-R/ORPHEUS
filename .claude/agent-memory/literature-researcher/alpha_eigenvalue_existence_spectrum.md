---
name: alpha-eigenvalue-existence-spectrum
description: Time (α) eigenvalue spectrum, existence of the fundamental α, Corngold limit, multigroup/SN artefacts, k-vs-α solving — extracted 2026-09-27; memo scratch/posing_sequence/alpha_existence/literature_memo.md
metadata:
  type: project
---

Extraction 2026-09-27 for the posing plan (#405). Full memo with page cites:
`scratch/posing_sequence/alpha_existence/literature_memo.md`.

Headlines (each `[M]` on the page unless marked):
- **Lehner–Wing (one-speed SLAB) proved a discrete α₀ ALWAYS exists** (Mockel 1966 NSE 26:279; Corngold 1964 NSE 19:85; B&G pp. 42/373; Carlvik 1968 Table I). ⛔ Dorning 2010 (Azmy–Sartori ch. 8, p. 406) mis-states it as "disappears for thin slab" — do not cite Dorning for that.
- **Bounded body, speeds bounded away from zero (one-speed, MULTIGROUP): NO continuous spectrum**; real eigenvalues exist below −min vΣ (Van Norton 1962; Ukai 1966; Carlvik 1968 NSE 31:295 Table II sphere to 14× vΣ, reproduced `[M]`; Betzler 2018 p.119; Larsen–Zweifel 1974 `[R]`). B&G p.374: "In the multigroup treatment, α₀ exists even for arbitrarily small systems".
- **Non-existence is a thermalization (v→0) phenomenon**: Albertoni–Montagnini 1966, Mockel 1966 (minimum thickness theorem, C_λ fictitious eigenvalue = the k(α) test), Corngold 1964 (max B*²), Vidav 1968 JMAA 22:144 (Thm III = Krein–Rutman fundamental).
- Physical kernels: lim_{v→0} vΣ = k₁ > 0 (Corngold Eq. 5); λ*=0 only in the constant-σ model, which has NO discrete modes (Corngold p.82).
- Discretised continuum: Ohanian–Daitch 1964 NSE 19:343 (energy-localised sign-changing modes); B&G p.375 "pseudo-fundamental"; McClarren 2019 p.862 α=−1.0216 artefact; Betzler 2018 excess eigenvalues reduce residual.
- k-vs-α: Hill 1983 (OSTI 6788588, α⁰=0 start, rebalance+Wielandt); B&G §1.5f + (9.8)/(9.13); Cullen 2003 UCRL-TR-201506 (α/v to RHS for α<0).
- ⚠ Mistral OCR was rate-limited (429); the 14 new PDFs have NO sidecars yet.
- NSE archive local; JMAA old papers downloadable via Elsevier API key in tools/research/elsevier.toml ([elsevier] api_key) with Accept: application/pdf.
