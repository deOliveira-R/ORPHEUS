---
name: gendf-ingest-truncation-and-n2n-probe
description: Durable facts about the GENDF ingest — EVERY channel is Legendre-truncated at ingest (elastic/inelastic/thermal at P2, (n,2n) at P0), load_isotope reads the .h5 cache so tape-only moments need _parse_gendf/_extract_mf6, the solver's min-over-materials L clamp, and the yield-strip hazard for ℓ≥1 — measured 2026-09-03 for the #426 probe.
metadata:
  type: project
---

Ground truth gathered for the #426 (n,2n)-anisotropy measurement (memo:
`scratch/_426_probe_ground.md`, untracked). The line numbers there are HEAD `1e02f6b1`;
the SHAPE below is what survives.

**The ingest truncates every channel, not just (n,2n).** `_build_isotope` keeps
`sig2_data[(0, 0)]` for MT=16 AND `_init_scattering`/`_accumulate_scattering` loop
`range(3)` for MT=2/51–91/221 — so `Isotope.sigS` is `[3][n_σ₀]` although the tape carries
NL=7 for elastic too (`cross_section_data.rst` "Three Legendre orders are always stored").
⟹ a "restore ℓ≥1 (n,2n)" experiment that leaves elastic at P2 mixes truncation orders
across channels; the honest ℓ=1,2-only arm needs no elastic re-read, the ℓ≤6 arm needs a
7-order re-read of MT=2 with σ₀ interpolation (elastic keys are `(ℓ, σ₀)`, (n,2n) keys are
`(ℓ, 0)` — σ₀-independent).

**`load_isotope` never touches the tape.** It reads `micro_xs/<name>.h5`, which stores only
what `Isotope` holds. Tape-only data comes from `_parse_gendf(_GXS_DIR/"X.GXS")` (~2 s,
275k cards for Be-9) + `_extract_mf6(mt, temp_idx, m)` → `(ifrom, ito, {(ℓ, iσ₀): vals})`,
1-based NJOY thermal-first indices; `tests/data/test_n2n_yield_convention.py` is the
worked example of driving these privately.

**Two silent traps for anyone augmenting `Mixture.SigS`:**
- `SNSolver` clamps `L = min(scattering_order, min(len(m.SigS)-1 over ALL materials))` —
  one 3-moment material drags a 7-moment probe down to L=2 with no message.
- `_strip_transfer_yield` is ℓ=0-only by contract (integer-yield admission on the row sums);
  ℓ≥1 matrices must be scaled by the ℓ=0 per-row factor `σ_MF3/rowsum_ℓ0`, never re-stripped.

**Convention (confirmed, complements [[harmonic-moment-field-units-convention]]):** the
`(2ℓ+1)` lives on the reconstruction `R` only (`LegendreBasis.reconstruct` /
`SphericalHarmonicBasis.reconstruct`); `Λ` applies bare `Σ_{s,ℓ}ᵀ`; GL-on-μ has Σw=2 so the
slab source is `½Σ(2ℓ+1)P_ℓΣ_{sℓ}ᵀφ_ℓ`. Elastic ℓ≥1 and (n,2n) ℓ≥1 are the same MF=6 record
type read by the same loop, so adding `2·Σ₂ₙ,ℓ` into `SigS[ℓ≥1]` is convention-consistent
by construction (the 2 is `N2NKernel.multiplicity`). ℓ≥1 never reaches
`SigT`/`absorption_xs`/`compute_keff`; the only off-sweep ℓ=1 reader is `transport_xs`
(diffusion/DSA).

**Feasibility `[M]`:** a 421-group, 20-cell, S8 reflected slab k-solve costs ~0.3 s/outer at
P0 and ~0.4 s/outer at P2 — condensation is NOT required for a 1-D probe.

**Two traps added 2026-09-03 (the #426 remedy blast-radius memo,
`scratch/_426_remedy_blast_radius.md`):**
- ⛔ **`_extract_mf6` PADS an NL=1 section to three keys with zeros** (`if n_lgn == 1:
  sig[(1,i)] += 0.0; sig[(2,i)] += 0.0`). A "distinct ℓ keys" census therefore reads **3**
  for NA023's MT=16 whose tape header says NL=**1** — I reported the pad as the tape on the
  first pass. Read the SEQ==1 header word (`m[i,2]`) or check `max|v|` of the ℓ≥1 keys before
  quoting an NL. `[M]` MT=16: NL=7 on 10 of 11 carrying tapes, NA023 NL=1 (padded), absent on
  B_010/H_001, **n_sig0 = 1 on all 11** (so a per-ℓ `sig2` is `[legendre]` only — no
  `interp_sig_s` analogue, unlike `sigS`'s `[legendre][sig0]`).
- The `*.h5` cache is **untracked and `.gitignore`d** (only the 13 `.GXS` are LFS-tracked);
  regeneration is the manual `orpheus/data/micro_xs/convert_gxs_to_hdf5.py` (README), no
  Makefile/test regenerates, and there is **no save/load round-trip test** — a schema change
  without a compatible loader makes recipes break on other checkouts while the h5-guarded
  tests SKIP (silent green). `load_isotope_h5` derives `n_legendre` from the `sigS/L{j}_S{k}`
  keys — the pattern a `sig2/L{j}` layout would follow.
- Dead-shim residue near the seam: `MaterialXSField.with_overridden_sig_s_and_n2n` has
  0 callers (its docstring claims `foldable_part`/`residual_part` — they build the field
  directly), and `derivations/diagnostics/diag_276_scattering_p0_fastpath_perf.py` calls
  `N2NMomentOperator(mat_xs=…, L=0)` — a signature that no longer exists.
