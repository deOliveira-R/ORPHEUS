# Literature Researcher — Memory Index

**§2 is a negative-lookup inventory** — "have I already extracted this
paper, and where?" It stays HOT because a miss costs a re-extraction.
Rows are telegraphic BY DESIGN: key + signature + file. Derivations and
value tables live in the linked file, greppable there; the row only
decides whether to open it.

## 1. Lessons

- [lessons.md](lessons.md) — READ FIRST. L1–L7, cited below as (L2), (L4).

## 2. Extraction inventory by topic

### Curvilinear / spherical SN
- [Slab & curvilinear LD](issue_158_linear_discontinuous_cell_update.md) — Larsen-Morel 1989 JCP 83 (4.1)-(4.3); diffusion-limit criterion LMM-1987 JCP 69 Table I; curvilinear LD UNPUBLISHED.
- [Multi-D LD closure](multi_d_ld_closure.md) — multi-D LD = BI/TRILINEAR DG-P1 (UBLD, MRM-2016 NSE), NOT simplex-P1 (Adams 2001: fails thick-diffusion limit).
- [Sphere-SN pole closure](sphere_sn_pole_closure_canonical.md) — Hébert §3.9.4 (3.418)-(3.439) IS canonical; `reduced_operator.py` Bailey docstring cites the WRONG paper (L3).
- [Carlson coupled pole](phase_d_carlson_coupled_pole.md) — Hébert (3.432)-(3.435) verbatim; µ=−1 seed on-the-fly; cyl analog (3.407)-(3.410).
- [Curvilinear sweep directness](curvilinear_sweep_directness_ruling.md) — RULING: DIRECT one-pass (Hébert / BMC / Lathrop-2000); NO published lagged-pole form; sole lag locus = albedo/white outer BC (3.415-3.417).
- [Central cell r=0 diamond](spherical_sn_central_cell_spatial_order.md) — Hébert §3.9.4 + Stacey §9.9 BOTH plain midpoint diamond; (3.430) φ = VOLUME-AVERAGE = the O(h) MMS gap.
- [Origin spatial order](sphere_sn_spatial_order_at_origin.md) — NO canonical ref gives an O(h²) central-cell SPATIAL closure (Lathrop 2000 + BMC = ANGULAR only); lead Wu 1999 NSE99-A2095.
- [Space ⊗ angle separability](space_angle_discretization_separability.md) — splits spatial (Larsen-Morel-Miller JCP 69 1987) vs angular (Bailey-Morel-Chang NSE 165 2010) → tensor-product architecture.

### SN solver / sweep theory
- [DSA primaries](dsa_primaries_full_extraction.md) — issue #2, 7 sources → `.claude/plans/dsa_literature_memo.md`: Alcouffe 1977 NSE 64:344 (print sign slips; ⅓ not ½) · Larsen I NSE 82:47 (four-step, K_N Table I) · McCoy-Larsen II NSE 82:64 · Morel 1982 NSE 82:34 · A&L 2002 · Adams-Martin 1992 M4S NSE 111:145 · Warsa 2004 k-IRAM. ⚠ NSE-147 TWIN trap; Σw=1 vs A&L Σw=2; absent WWM-147(2) / Yavuz-Larsen'88.
- [Adams & Larsen 2002](adams_larsen_2002_iteration_theory.md) — PNE 40 FULL → `.claude/plans/phase_i_survey_adams_larsen_2002.md`: ρ_SI=c (2.17/3.30), DSA 0.2247c (2.50/3.65), Krylov Tables 1/2; traps (λ in mfp; cyl µ=angle; OCR ⅓→½).
- [Morel 1989 SI ≡ apply](morel_1989_si_vs_apply_equivalence.md) — p.75: SI sweep ≡ the linear system IFF angular leakage is lower-triangular; else Jacobi-on-RHS, ρ→1.
- [Larsen & Morel 2010 review](larsen_morel_2010_review_extraction.md) — 84 pp; NO cyl-SN / adjoint / k-posing / verification; ρ=c + DSA 0.23c; their β ≡ BMC α, α = SPATIAL WD; PDF offset +12.
- [Method naming evidence](method_canonical_naming_evidence.md) — quotes pinning "singular eigenfunction expansion" (Mitsis / WM-72 / Atalay / Sanchez-1977) vs "trajectory-based/resolvent" (Sanchez 2002/1986, PS-1982).

### Multi-group source indexing
- [χ fission-source indexing](peierls_mg_fission_source_chi_indexing.md) — Hébert (3.57)/(3.58): emission χ is purely LOCAL — shares the SOURCE point with νΣ_f and φ.

### Energy condensation / group collapse
- [Condensation formulas](energy_condensation_collapse_formulas.md) — Hébert §3.5 (3.103)-(3.105)/(3.112 χ) ≡ Stammler Ch.VI(6a-6d); MALOCS REQUIRES nesting, GROUPR/OpenMC re-integrate continuum; NJOY IWT; GEC (Rahnema-Douglass-Forget 2008 NSE 160:41) = projection, rank-0 = flux-wt avg. Memo `.claude/plans/p5_condensation_literature.md`; #274.
- [Adjoint bilinear collapse](adjoint_bilinear_collapse_p6.md) — P6 #281. B&G Ch.6 [LOCAL]: (6.135) ⟨ψ†σψ⟩; (6.136) scattering = per-pair SINK-adjoint × SOURCE-flux; carriers (6.125)/(6.126); fission dyad χ†_g·(νσ_f)_{g′}; (6.90), δk (6.71). + Hébert Ch3 / Dorning ch8 / Roy ch4. Memo `.claude/plans/p6_literature_memo.md`.

### SH angular flux → moment projection
- [SH projection root cause](sh_flux_moment_projection_root_cause.md) — "M exists BECAUSE of anisotropic scattering" STANDS. Hébert (3.55)→(3.54), Brockmann (47), Ahrens LDO (7), Fletcher (7) → HarmonicFrame on Scattering, NOT phase-space.

### Interface currents (rank-N / DP_N)
- [Sanchez-McCormick 1982](rank_n_interface_current_canonical.md) — §III.F = the LONE rank-N source in the corpus; Ligou / Sanchez-2002 / Stamm'ler / Stacey all scalar DP-0.
- [S-M 1982 Eq.166](sanchez_mccormick_rank_n_per_face.md) — Eq.166 cluster verbatim; the σ_t-dependent N=2 residual = a missing (Ω·n) µ-weight, NOT (ρ/R)².
- [Rank-N 4-ref synthesis](rank_n_closure_four_references_synthesis.md) — Ligou 1982 + Sanchez 2002 + Stamm'ler Ch.IV + Stacey 2007 all scalar/DP-0; ORPHEUS F.4 IS the textbook closure.
- [Hébert 2009 Ch.3 IC](hebert_2009_ch3_interface_currents.md) — DP_N matches §III.F.1; only 2-D Cartesian instantiated; 1-D cyl/sph → scalar F.4; Villarino-Stamm'ler renorm.
- [Stamm'ler 1983 Ch.VI](stammler_1983_ch6_interface_currents.md) — Ch.VI is SN; Ch.IV is CP + scalar (rank-0) IC. NO rank-N IC anywhere in Stamm'ler-Abbate.
- [Sanchez 1977 NSE 64](sanchez_1977_nse64_canonical_ic.md) — PDF read. Rank-3 = 3 modes/face; NO multi-collision K_inf (fixed-source only); Class B cap ≈3% (Tables VI-VII).
- [DP_N curvilinear gap](dpn_curvilinear_nonexistence.md) — no Stepanek-style DP_N k_eff tables for cyl/sphere; Sood 2003 F_N is the correct truth set.
- [Rank-N curvilinear leads](rank_n_ic_curvilinear_literature_leads.md) — Bogado Leite 1998, Krishnani 1982, Mohanakrishnan 1982 = the only ANE title matches; Bogado Leite orphaned — read before any novelty claim.
- [Lambert ↔ Marshak](direction_q_lambert_marshak_derivation.md) — Sanchez 2014 NSE 177 is gauge-theoretic at differential P_N, NOT integral CP; solid-harmonic change-of-basis on [0,1].
- [Knyazev 1993](knyazev_1993_cylinder_anisotropic_ic.md) — *Atomic Energy* 74(5):368-374, DOI 10.1007/BF00844623 = THE cylinder Ki_{2+k} 3-D source; "Knyazev-Selivanov 2014" = phantom (L2).

### Collision probability (slab/cyl/sphere kernels)
- [CP moment integrals](cp_moment_integrals.md) — J_k slab + T_k cylinder closed forms: A&S §5.1 + Hébert §3.4-3.5. Stamm'ler Ch.IV is rank-0 flat-source only.
- [Carlvik 1967](carlvik_1967_finite_cylinder_cp.md) — PDF read. 2-pp Note, E_3 kernels (cuboid + cylinder); NO Bickley / rank-N. Real Wigner-Seitz CP paper = Carlvik 1965 Geneva Vol.2 p.225.
- [Sphere CP + white BC](bickley_naylor_sphere_white_bc.md) — sphere kernel = E_2, NOT Bickley/Ki_n (cylinder-specific); Hébert §3.8.5 (3.324)-(3.336); G_bc = 4·P_esc.
- [Bickley libraries](bickley_function_libraries.md) — Davierwalla 1982 / Lorensi 2025 / Amos 1983 TOMS 609 vs ORPHEUS `ki_n_float` (78 µs; 1e-13). Atkinson F_N floor uses E_1.
- [Cylinder Peierls](phase4_cylinder_peierls.md) — Ki1 vs Ki3, 1/π prefactor, chord branches, rank-N_y white BC, benchmark radii.

### Continuous-µ specular (Peierls Green's function)
- [Continuous-µ specular](phase5_continuous_mu_specular.md) — Hébert §3.8.3 + Sanchez 2002 Eq.15 = the only textbook multi-bounce forms; Sanchez 1986 TTSP / PS-1982 / Milgram 1978 = analytical refs.
- [Sanchez 1986 sphere](phase5_sanchez_1986_sphere_specular.md) — PDF read. Eq.(A6) IS the continuous-µ multi-bounce specular kernel; T(µ)=1/(1−α·exp(−2aµ)); multi-region via piecewise τ(µ).
- [Peierls Green's function](peierls_greens_function_lit.md) — Sanchez 2002 §III is LATTICES not sphere; Sanchez 1986 supersedes PS-1982; use Variant α (closure-free op on A6); no sphere-specular k_eff precedent.
- [Peierls extensions wishlist](peierls_greens_extensions_lit.md) — expansion list across geom/topology/BC/MG/multi-region; 22 PDFs LOCAL; wishlist Garcia 2021 JCP + Knyazev 1993 + Hébert Ch.6-9.
- [PS-1982 + Garcia 2021](ps1982_and_garcia_extraction.md) — PS-1982 Eq.(21) = structurally-indep vacuum-sphere kernel (NO k_eff → homog-eigval + Sood c_crit); Garcia 2021 fixed-source ONLY, Tables 5/12/18.
- ["Green's function" 3 meanings](sanchez_chandrasekhar_three_meanings.md) — M1 Chandrasekhar angular G / M2 Sanchez source-reduction kernel / M3 trajectory MoC; cylinder M2 needs Carlvik 1965 / Knyazev 1993.

### F_N method (criticality benchmarks)
- [Sood LA-13511 1999](sood_fn_method_full_extraction.md) — 75-problem catalogue + naming convention; F_N machinery is in KLL 1974 / Westfall-Metcalf / Siewert-Thomas, NOT the report.
- [Sood 2003 vs 1999](sood_2003_vs_1999_extraction.md) — the journal port adds NO method beyond a general-multigroup k_inf appendix; both are TEST SET, not method papers (L4).
- [Sood 2003 cylinder values](sood_2003_cylinder_benchmarks.md) — `Ua-1-O-CY` (c=1.30) r_c=1.72500292 mfp; `PUb-1-O-CY` (c=1.40) r_c=1.396979 mfp; Westfall F_N, vacuum-only.
- [Kaper-Lindeman-Leaf 1974](kaper_lindeman_leaf_1974_fn_method.md) — KLL NSE 54-94 paywalled, NO OA mirror; derivation also in Garcia-Siewert 1979 (paywalled); pivot to Case singular-EF (B&G Ch.2.6).
- [Siewert-Thomas 1986 sphere](sphere_fn_method_extraction.md) — sphere = slab with anti-symmetric BC (one sign flip); F_α recursion + Chebyshev collocation + X-function reusable from slab.
- [F_N flux reconstruction](meaning1_already_in_fn_method.md) — KLL X-function + A(ν)/B(ν) Fredholm + `*_angular_flux_from_scalar` ALREADY present (slab+sphere isotropic); gap = API exposure. See `fn_method/{slab,sphere}/flux_reconstruction.py`.
- [Atkinson product-Nyström](atkinson_product_nystrom_log_kernels.md) — 1997 book §4.2 + 1972 Numer. Math. (both LOCAL); Bickley→log via Ki1(τ)~τ log τ; product-Simpson + q=4 graded mesh = the F_N fix.

### Anisotropic / reflected criticality (Case-EF & Galerkin)
- [Atalay 1997](atalay_1997_reflected_anisotropic.md) — PNE 31(3):229-252, NOT in Zotero. Case singular-EF + Fredholm, unified slab/sphere; Tables 2-10 = reflected+anisotropic benchmark suite. NOT F_N (L5) → `case_method/`.
- [Burkart-Ishiguro-Siewert 1976](burkart_ishiguro_siewert_1976_two_region_anisotropic.md) — NSE 61:72-77, PDF read. Case-EF + Chandrasekhar H; Tables I (Milne endpoints) + II (critical half-thicknesses), 40 rows each. NOT F_N.
- [Dahl-Sjöstrand 1979](dahl_sjostrand_1979_anisotropic_slab_sphere.md) — NSE 69:114-125. Legendre-Galerkin on Carlvik's integral eq; full c-eigenvalue SPECTRUM (≤11) × {µ̄, d}; PRIMARY for P_1 slab+sphere. NOT F_N.
- [Neshat-Maiorino 1980](neshat_maiorino_1980_reflected_slab_fn.md) — ANE 7:79-81, PDF read. Reflected-slab F_N (3 boundary integrals → 3(N+1)); unlocks ~10-14 Sood reflected cases; Table 1 τ_c.

### Corpus bibliography (#231 Phase G2)
- [docs/refs.bib ledger](refs_bib_g2_corrections.md) — COMPLETE, 59 entries, pybtex-clean, DOIs CrossRef-verified. **Trust the BIB over any page/docstring definition**; keys immutable, only FIELDS fixed (per-entry ledger in file). Swaps: PS1982→PomraningSiewert1982, Sood1999→SoodLA13511_1999.

### Tooling
- [Zotero failover](reference_zotero_flakiness.md) — dead-server signature = 0 hits on known items + conn-refused on port 23119 → Tier 2 (L7).

## 3. Active state

- **No open literature pull.** Only open SN issue with literature relevance: **#236** (spatial ⊗ angular product) — basis [[space_angle_discretization_separability]].
- Campaign mentions in topic files (#158/#240, #257/ERR-063, #196) are CLOSED archaeology kept for their equations — reconcile any "in-flight" claim against git, never a frozen note.
- "Zotero DOWN this session" notes in topic files are EPHEMERAL — re-query when an item's annotations matter.
</content>
