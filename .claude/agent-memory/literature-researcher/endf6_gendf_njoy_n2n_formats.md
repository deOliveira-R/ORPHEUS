---
name: endf6-gendf-njoy-n2n-formats
description: The ENDF-6 + NJOY/GENDF nuclear-data FORMAT specs — acquisition URLs, page offsets, and the settled (n,2n)/MT=16 answers (LAW menu, GENDF record layout, Legendre convention, IG=0 fission-only, Be-9 = LAW=7 LAB)
metadata:
  type: reference
---

# ENDF-6 / GENDF / NJOY — the nuclear-DATA-FORMAT corpus

⭐ **This is a whole domain the §2 inventory previously had NO entry for.** Prior to
2026-08-31, `scratch/literature/` (78 PDFs) contained **zero** nuclear-data-format
documents — it is a *transport-theory* library. Both primary specs are **free, born-digital
(real text layer, not scans), and fetch in seconds from fixed URLs.** Do not treat a
data-format question as a Tier-2 literature hunt.

## Acquisition — both specs, one curl each (⭐ verified working 2026-08-31)

```bash
# ENDF-6 Formats Manual — CSEWG ENDF-102 / BNL-203218-2018-INRE, rev 215, 418 pp
curl -sSL -o ENDF-102-2018.pdf https://www.nndc.bnl.gov/endfdocs/ENDF-102-2018.pdf
# NJOY2016 manual — LA-UR-17-20093, upd. NJOY2016.53, 816 pp  (raw, NOT the releases API)
curl -sSL -o njoy16.pdf https://raw.githubusercontent.com/njoy/NJOY2016-manual/master/njoy16.pdf
# any ENDF/B-VIII.0 evaluation, by MAT:
curl -sSL -o x.zip https://www-nds.iaea.org/public/download-endf/ENDF-B-VIII.0/n/n_0425_4-Be-9.zip
```
⛔ The GitHub **releases** URL 404s; use `raw.githubusercontent.com`. The IAEA `.dat`
direct path 404s; only the `.zip` works. MAT is zero-padded 4 digits and the element
string must match exactly (`n_0528_5-B-11`, **not** `n_0525_5-B-11` — B-10 is 525).

## Page offsets (established 2026-08-31; both are stable published revisions)

| doc | mapping |
|---|---|
| ENDF-102 (2018) | **printed p. = PDF p. − 18** |
| NJOY2016 (LA-UR-17-20093) | **printed p. = PDF p. − 14** |

## ⛔⛔ The search hazard that WILL bite you again

`grep` here is **ugrep with `-I`**, which classifies a `pypdf`-extracted `.txt` as
**binary and silently skips it — 0 hits, no error, indistinguishable from a clean file.**
I lost a cycle to this. **Extract to text with page markers and search in PYTHON**:

```python
parts = re.split(r'@@@PDFPAGE (\d+)@@@', txt)
pages = {int(parts[i]): parts[i+1] for i in range(1, len(parts), 2)}
```
(Helper left at `scratch/specs/s.py`.) Same family as the L-family filter lessons — a
broken filter and an empty corpus print the same thing. See [[lessons]] L-015.

## ⚠ pypdf breaks FRACTIONS across lines — verify any normalization visually

`(2ℓ+1)/4π` extracted as `2ℓ + 1 / 4π` on two lines. The `/2` vs `/4π` distinction is
exactly what a Legendre-convention question turns on. **Rendered-page check is mandatory
for any normalization constant**, even though these are born-digital PDFs where the
ERR-032 OCR hazard does not apply. Both were verified 2026-08-31 (L-010's discipline,
new mechanism: layout loss, not OCR error).

---

# Settled answers — MT=16 / (n,2n). Full memo: `scratch/n2n_data_format_spec.md` (855 lines)

## The Legendre convention (⭐ the most reusable fact)

- **ENDF-102 Eq. (6.3)**, printed p. 135: `f = Σ_{l=0}^{NA} [(2l+1)/2] f_l P_l`, with the
  printed note that `f_0` is the **angle-integrated total** ≡ File 5's `g(E,E′)`.
- **NJOY Eq. (242)**, printed p. 188: `σ_X = Σ_ℓ [(2ℓ+1)/4π] σ_Xℓ P_ℓ(µ₀)`.
- The `/2` vs `/4π` is a **units difference, not a conflict**: ENDF's `f` is per unit
  cosine and enters via Eq. (6.1) `σ_i = σ·y·f/2π` [b/sr]. Combining gives
  **`σ_Xℓ = σ(E)·y(E)·f_ℓ`**, hence two consequences worth memorising:
  1. **`Σ_g′ σ_{ℓ=0,g→g′} = σ(E)·y(E)`** — the MULTIPLICITY IS ALREADY FOLDED IN.
     `[EVAL]` measured **2.000000** on Be-9 MT=16 over 50 groups. Never multiply by 2 again.
  2. **`σ₁/σ₀ = µ̄` identically** — yield and cross section cancel. No `(2ℓ+1)` rescale
     on read. (Elastic control: `→ 2/(3A)`.)

## GENDF record layout (NJOY §8.17, printed pp. 228–229)

- header MF=1/MT=451: CONT `ZA, AWR, 0, NZ, −1, NTW` (the literal **−1** = "this is GENDF")
  + LIST `TEMPIN, 0, NGN, NGG, NW, 0` / `TITLE, SIGZ, EGN, EGG`, `NW = NTW+NZ+NGN+1+NGG+1`.
  FEND follows; **SEND is omitted**. MF/MT ordering requirements relaxed.
- reaction CONT: `ZA, AWR, NL, NZ, LRFLAG, NGN`
- per-group LIST: `TEMP, 0, NG2, IG2LO, NW, IG`, **`NW = NL·NZ·NG2`**; last record is `IG=NGN`,
  given even if zero.
- transfer matrix `A` = **`FLUX(NL,NZ), MATRIX(NL,NZ,NG2−1)`** — the leading `NL·NZ` "flux
  words" are the Eq. (251) weighting denominator, not padding; ℓ varies fastest.
  In Bondarenko mode all `NL` flux entries are **identical** (printed p. 256).

## The five settled MT=16 rulings

1. **Home**: MF=6 recommended, MF=4+MF=5 still legal (ENDF §5.1 names MT=16 explicitly).
   `LAW=1` is the (n,2n) law (multi-body); `LANG=1` Legendre / `LANG=2` Kalbach–Mann /
   `LANG=11-15` tabulated. **`NA` is NOT an anisotropy order until LANG is read** — it is
   the Legendre order only under LANG=1; under LANG=2 it is a count ≤2 of Kalbach params.
2. ⛔ **`LAW=6` (phase space) is EXACTLY isotropic in CM** (Eq. 6.21 has no µ) but its LAB
   form (Eq. 6.25) is µ-dependent. **CM-isotropic ⇏ LAB-isotropic.**
3. ⛔⛔ **The processing-artifact trap: NJOY Eq. (361)** builds LAB moments as
   `p_Lℓ = ∫ p_C · P_ℓ(µ_L) · J dµ_L` — `P_ℓ` at the LAB cosine, `p_C` at the CM cosine, so
   **ℓ≥1 moments are non-zero even for an isotropic CM distribution.** ⟹ for any `LCT=2`
   evaluation a non-zero `σ₁/σ₀` in GENDF **mixes evaluated anisotropy with kinematic
   boost and cannot be decomposed from the GENDF alone.** Only `LCT=1` files are clean.
4. **`IG=0` is FISSION-ONLY** (NJOY printed p. 229): emitted only for the fission matrix and
   fission/capture photon-production matrices, via the §8.6 break-energy machinery
   (`σ_f,g′→g = σ^HE + χ^LE_g·σ^LE_fPg′`, Eq. 283). ⟹ **MT=16 is structurally incapable of
   carrying one.** The separable χ is a **GROUPR construct, not an ENDF object** — GROUPR
   *detects* separability by scanning File 5 for a break energy, and looks nowhere but fission.
5. **No ENDF representation stores an E-INDEPENDENT secondary spectrum** — every File 5 law
   is E-dependent (LF=1 TAB2 over E; LF=5/7/9/11/12 via `θ(E)`, `a(E)`, `b(E)`). Nearest is
   LF=5's fixed `g(x=E′/θ(E))`, which NJOY restricts to delayed neutrons.

## ⭐ Per-nuclide census (ENDF/B-VIII.0), 2026-08-31 — reuse before re-fetching

| nuclide | MT=16 | LCT | angular content |
|---|---|---|---|
| **Be-9** (MAT 425) | MF=6 **LAW=7** | **1 = LAB** | 24 E × **21 µ** × ≤93 E′; NK=2 (n, **⁴He** y=2); **94 % of the whole file**; `[EVAL]` µ̄ = +0.58 at threshold → **+0.25…+0.27** plateau, f(+1)/f(−1) = 4–12; ⟨E′⟩ 0.051→6.216 MeV (**×122** ⟹ rank-1 impossible) |
| Na-23 (1125) | **MF=4 + MF=5** | 1 | ⭐ MF=4 is **2 lines: LTT=0, LI=1 = ALL ISOTROPIC**; MF=5 LF=9 evaporation, U=12.414 MeV |
| B-11 (528), O-16 (825) | MF=6 LAW=1 **LANG=2** | 2 = CM | Kalbach–Mann, NA=1 (`f_0`+`r` only — higher ℓ are ANALYTIC expansion, not evaluated) |
| Zr-90 (4025) | MF=6 LAW=1 LANG=1 | 2 = CM | Legendre, **NA→14** (richer than lord=6 keeps) |
| U-235 (9228) / U-238 (9237) | MF=6 LAW=1 LANG=1 | 2 = CM | Legendre, NA→4 / **NA→15**; NK=3 |
| H-1 (125), B-10 (525) | **absent** | — | no MT=16 section at all |

## What controls GENDF `NL` — DATA-DEPENDENT, and the manual alone cannot tell you

`lord` (GROUPR card 2) is a **ceiling**; `subroutine init` picks a per-`(mfd,mtd)` default
(printed p. 255, **table not printed**); `getco` truncates File-4 data at `toler=1e-6`.
⟹ **`NL=7` DOES prove the evaluation is not declared isotropic** — proof by the two
counterexamples under the same `lord=6` run: **Na-23 → NL=1** (evaluation says `LI=1`) and
**MT=18 → NL=1** (GROUPR's OWN assumption — §8.6 opens *"Assuming isotropy"*).
⚠ But `NL=7` does **not** prove 7 *evaluated* moments (LANG=2 nuclides supply only `f_0`+`r`;
Zr/U are being truncated), and within a section **all NL moments are written including zeros**.
Also: MF=6/MT=16 has **`NZ=1`** (not self-shielded) while MF=6/MT=2 has `NZ=6`.

## Production practice (Q5) — splits by LIBRARY FORMAT, no single standard

- **RETAIN per-ℓ**: ISOTXS/CCCC — (n,2n) is a first-class block, `IDSCT=300+NN`, own `LORDN`;
  *"total scattering matrix is the sum of the elastic, inelastic, and (n,2n) matrices"*
  (NJOY printed p. 423). `SN2N` = **reaction** xs; the matrix sums to **2×** = production xs
  (printed p. 420) — same convention as GENDF.
- **COLLAPSE to P0**: the **WIMS library format** — *"non-thermal **P0 scattering matrix** …
  summing over all of the reactions … transport corrected by subtracting the … P1 matrix"*
  (NJOY printed p. 562). And the CCCC **diffusion** transport xs, Eq. (469)
  `σ_tr = σ_t1 − Σσ_e1`, subtracts **ELASTIC P1 ONLY** ⟹ (n,2n) anisotropy dropped.
- ⭐ **The common SILENT approximation**: Boyd, Nelson, Romano, Shaner, Forget, Smith (2019),
  *Nucl. Technol.*, **DOI 10.1080/00295450.2019.1571828** (CrossRef-verified, 52 cites),
  §III.D–E: the scattering-multiplicity matrix `υ̂` is built from **ℓ=0 only** (Eq. 13) and
  then multiplies **every** ℓ as a scalar (Eq. 18) ⟹ **a P0 treatment of (n,2n) SHAPE wearing
  a P_N matrix.** GENDF is strictly richer. Serpent 2 (`serpent.vtt.fi/docs/user_guide/
  gc_generation.html` §2.2.5) gives `Σ_sp = Σ_s + 2Σ_2n + 3Σ_3n + …` and P_n to n=7 but
  **explicitly disclaims n>0** as "not fully consistent with deterministic transport theory".
- ⛔ **No primary source anywhere argues (n,2n) anisotropy is negligible on MAGNITUDE
  grounds.** P0 is inherited from the whole matrix being P0, never justified for MT=16
  specifically. ENDF §6.3.3 argues the *opposite* for evaluators.

## Open / unclosed

1. ENDF §6.3.3's identical-particle CM-symmetry clause (printed p. 150) — scope UNRESOLVED:
   does it bind every multi-body exit channel with two like particles, or only
   identical-particle *elastic* scattering? No disambiguating statement in ENDF-102.
2. The `init` (mfd,mtd)→NL table needs the GROUPR **source**, not the manual.
3. CASMO / APOLLO2-3 / DRAGON MT=16 Legendre treatment — no primary doc obtained (try HAL
   for APOLLO/DRAGON; CASMO is proprietary).
4. JEFF-3.3 / JENDL-5 Be-9 MT=16 not checked.

Artifacts: `scratch/specs/{ENDF-102-2018.pdf,endf102.txt,njoy16.pdf,njoy16.txt,n_*.dat,
openmc_mgxs.pdf,serpent_gc.html,s.py}`. Related: [[lessons]] L-015/L-016.
