# Literature Researcher — Memory Index

Slim index. Section 1 = behavioral lessons (read first). Section 2 is
the BULK: a "what have I already extracted, and where" lookup over
durable paper-extraction RECORDS (equations + numbers, notation maps,
method-classification verdicts, reference-value tables). These stay
useful indefinitely — a paper's extracted equations do not expire.
Section 3 = git-true active state.

## 1. Lessons (read first)

- [lessons.md](lessons.md) — 7 research-discipline lessons: Bailey µ
  geometry trap (L1); phantom citation = memory-only provenance (L2);
  verify a docstring's citation (L3); catalogue ≠ method source (L4);
  classify method by reading the paper not its citing context (L5);
  local-folder-first then Zotero then Tier-2 (L6); dead-Zotero failover (L7).

## 2. Durable reference — extraction inventory by topic

### Curvilinear / spherical SN (discretization, pole closure)
- [issue_158_linear_discontinuous_cell_update.md](issue_158_linear_discontinuous_cell_update.md) — slab LD = Larsen-Morel 1989 JCP 83 Eqs (4.1)-(4.3); diffusion-limit criterion = LMM-1987 JCP 69 Table I; curvilinear LD UNPUBLISHED (derive by fusing slab-LD + M-M curv-DD).
- [multi_d_ld_closure.md](multi_d_ld_closure.md) — Cartesian multi-D LD = BILINEAR/TRILINEAR DG-P1 (UBLD, 2^d corner moments, MRM-2016 NSE), NOT simplex-P1; Adams 2001 proved simplex fails thick-diffusion-limit on quads.
- [sphere_sn_pole_closure_canonical.md](sphere_sn_pole_closure_canonical.md) — Hébert §3.9.4 (3.418)-(3.439) IS the canonical sphere-SN pole stencil; the `reduced_operator.py` Bailey docstring cites the WRONG paper (see L3).
- [phase_d_carlson_coupled_pole.md](phase_d_carlson_coupled_pole.md) — Hébert (3.432)-(3.435) coupled-pole sweep verbatim; inward µ=−1 seed is on-the-fly; outward cascade reads full radial profile. Cyl analog (3.407)-(3.410).
- [curvilinear_sweep_directness_ruling.md](curvilinear_sweep_directness_ruling.md) — RULING: 1-D curv-SN sweep is DIRECT one-pass (Hébert/BMC-17,42,55/Lathrop-2000 "sequential"); NO published lagged-pole formulation (lag = verdict (c)); only std lag locus = albedo/white OUTER BC (Hébert 3.415-3.417 shooting = direct alt); α-recursion ×2/sign normalization table.
- [spherical_sn_central_cell_spatial_order.md](spherical_sn_central_cell_spatial_order.md) — Hébert §3.9.4 + Stacey §9.9 BOTH use plain arithmetic-midpoint diamond at r=0; Hébert (3.430) defines φ as VOLUME-AVERAGE (the O(h) MMS gap source).
- [sphere_sn_spatial_order_at_origin.md](sphere_sn_spatial_order_at_origin.md) — NO canonical ref gives an O(h²) central-cell SPATIAL closure; Lathrop 2000 + Bailey-Morel-Chang study ANGULAR error only. Remedy lead = Wu 1999 NSE99-A2095.
- [space_angle_discretization_separability.md](space_angle_discretization_separability.md) — diffusion-limit lit SPLITS into spatial (Larsen-Morel-Miller JCP 69 1987) + angular (Bailey-Morel-Chang NSE 165 2010) per-axis conditions; factorize but jointly required → tensor-product architecture.

### SN solver / sweep theory
- [adams_larsen_2002_iteration_theory.md](adams_larsen_2002_iteration_theory.md) — A&L PNE 40 FULLY extracted → `.claude/plans/phase_i_survey_adams_larsen_2002.md`: ρ_SI=c (2.17/3.30), DSA 0.2247c (2.50/3.65), consistency+four-step, Krylov Tables 1/2 + p.140 GMRES-outside-DSA quote, 14 M2 triples; traps (λ in mfp units; cyl µ=angle; OCR ⅓→½; SPI [131] slip).
- [morel_1989_si_vs_apply_equivalence.md](morel_1989_si_vs_apply_equivalence.md) — Morel 1989 p.75: SI sweep ≡ discrete-ordinates linear system iff angular leakage matrix is lower-triangular; full matrix → Jacobi-on-RHS, ρ→1.
- [larsen_morel_2010_review_extraction.md](larsen_morel_2010_review_extraction.md) — L&M-2010 84-pp review surveyed for the doc corpus: NO cyl-SN/adjoint/k-posing/verification; ρ=c + DSA 0.23c + thick-limit eq map; their β≡BMC α, their α=SPATIAL WD; book-PDF offset +12.
- [method_canonical_naming_evidence.md](method_canonical_naming_evidence.md) — verbatim quotes pinning "singular eigenfunction expansion" (Mitsis/WM-72/Atalay/Sanchez-1977) and "trajectory-based/resolvent" (Sanchez 2002/1986, PS-1982).

### Multi-group source indexing
- [peierls_mg_fission_source_chi_indexing.md](peierls_mg_fission_source_chi_indexing.md) — Hébert (3.57)/(3.58): fission emission χ is purely LOCAL, shares the SOURCE point with νΣ_f and φ; χ is a property of the fissioning nuclide.

### Energy condensation / group collapse (vectors, 2-axis scattering, χ; nesting; IWT; projection)
- [energy_condensation_collapse_formulas.md](energy_condensation_collapse_formulas.md) — Hébert §3.5 (3.103)/(3.104 src-flux-wt sink-summed)/(3.105)/(3.112 χ plain-sum) ≡ Stammler Ch.VI(6a-6d) [both LOCAL]; MALOCS REQUIRES nesting (correspondence array) vs GROUPR/OpenMC re-integrate continuum (any structure); NJOY IWT w(E) taxonomy; GEC (Rahnema-Douglass-Forget 2008 NSE160:41)=projection whose rank-0=flux-wt avg. Memo `.claude/plans/p5_condensation_literature.md`. For #274.
- [adjoint_bilinear_collapse_p6.md](adjoint_bilinear_collapse_p6.md) — P6 #281 COMPLETE: B&G Ch6 [LOCAL] = the per-channel bilinear source — (6.135) vector ⟨ψ†σψ⟩ per-moment; **(6.136) scattering = per-pair SINK-adjoint×SOURCE-flux**; carriers (6.125)/(6.126): plain-integral flux + FLUX-WEIGHTED-avg adjoint ⟹ ⟨φ*φ⟩_G=Φ*_G·Φ_G ((a)≡(b) under separability); fission dyad factored χ†_g·(νσ_f)_{g′}; theorem (6.90) + δk (6.71) + §6.2c consistency. Plus Hébert Ch3 / Dorning ch8 / Roy ch4 partition. Memo `.claude/plans/p6_literature_memo.md`.

### SH angular flux→moment projection (M/R) — root cause
- [sh_flux_moment_projection_root_cause.md](sh_flux_moment_projection_root_cause.md) — claim "M (ψ→φ_ℓ) exists BECAUSE of anisotropic scattering" STANDS (no falsifier). Hébert (3.55)→(3.54), Brockmann (47), Ahrens LDO (7)+abstract; PN basis DIAGONALIZES scattering (Fletcher (7)), streaming tridiagonal. FCS Eq.(6) q_1st=HΨ_u is scattering op; external/BC sources SPECIFIED in SH (input) or ordinate-space. → HarmonicFrame on Scattering operator, NOT angular phase-space.

### Interface currents (rank-N / DP_N), curvilinear
- [rank_n_interface_current_canonical.md](rank_n_interface_current_canonical.md) — Sanchez-McCormick 1982 §III.F is the LONE rank-N source in the textbook corpus; Ligou/Sanchez-2002/Stamm'ler/Stacey all use scalar DP-0.
- [sanchez_mccormick_rank_n_per_face.md](sanchez_mccormick_rank_n_per_face.md) — Eq.166 cluster verbatim; the σ_t-dependent N=2 residual = a missing (Ω·n) µ-weight, not (ρ/R)².
- [rank_n_closure_four_references_synthesis.md](rank_n_closure_four_references_synthesis.md) — Ligou 1982 + Sanchez 2002 + Stamm'ler 1983 Ch.IV + Stacey 2007 all scalar/DP-0; ORPHEUS F.4 IS the textbook closure; plateau = non-conservation of higher modes.
- [hebert_2009_ch3_interface_currents.md](hebert_2009_ch3_interface_currents.md) — DP_N machinery matches §III.F.1; instantiated only 2D Cartesian; 1D cyl/sph reduces to scalar F.4; admits conservation failure → Villarino-Stamm'ler renorm.
- [stammler_1983_ch6_interface_currents.md](stammler_1983_ch6_interface_currents.md) — Ch.VI is SN; Ch.IV is CP + scalar (rank-0) IC. Stamm'ler-Abbate has NO rank-N IC anywhere.
- [sanchez_1977_nse64_canonical_ic.md](sanchez_1977_nse64_canonical_ic.md) — PDF read. Rank-3 = 3 modes/face; NO multi-collision K_inf (fixed-source only); rank-3 cap on heterogeneous Class B ~3% (Tables VI-VII).
- [dpn_curvilinear_nonexistence.md](dpn_curvilinear_nonexistence.md) — no Stepanek-style DP_N k_eff tables for cyl/sphere; Sood 2003 F_N is the correct truth set.
- [rank_n_ic_curvilinear_literature_leads.md](rank_n_ic_curvilinear_literature_leads.md) — Bogado Leite 1998, Krishnani 1982, Mohanakrishnan 1982 = only ANE title matches for rank-N IC in cyl/sphere; Bogado Leite orphaned (read PDF before novelty claim).
- [direction_q_lambert_marshak_derivation.md](direction_q_lambert_marshak_derivation.md) — Sanchez 2014 NSE 177 is gauge-theoretic at differential P_N (not integral CP); solid-harmonic Lambert↔Marshak change-of-basis on [0,1] is the best symbolic experiment without PDFs.
- [knyazev_1993_cylinder_anisotropic_ic.md](knyazev_1993_cylinder_anisotropic_ic.md) — A.P. Knyazev 1993 *Atomic Energy* 74(5):368-374, DOI 10.1007/BF00844623 — THE cylinder Ki_{2+k} 3-D source; "Knyazev-Selivanov 2014" was a phantom (see L2).

### Collision probability (slab/cyl/sphere kernels, moments)
- [cp_moment_integrals.md](cp_moment_integrals.md) — J_k slab + T_k cylinder closed forms: A&S §5.1 + Hébert §3.4-3.5. Stamm'ler Ch.IV is rank-0 flat-source only.
- [carlvik_1967_finite_cylinder_cp.md](carlvik_1967_finite_cylinder_cp.md) — PDF read. 2-pp Note, E_3 kernels for finite cuboid + cylinder; NO Bickley-Naylor / boundary closure / rank-N. The real Wigner-Seitz CP paper = Carlvik 1965 Geneva Vol.2 p.225.
- [bickley_naylor_sphere_white_bc.md](bickley_naylor_sphere_white_bc.md) — sphere CP canonical kernel = E_2 (NOT Bickley/Ki_n, which are cylinder-specific); Hébert §3.8.5 (3.324)-(3.336); G_bc(r)=4·P_esc(r) by reciprocity.
- [bickley_function_libraries.md](bickley_function_libraries.md) — Davierwalla 1982 vs Lorensi 2025 vs Amos 1983 TOMS 609 vs ORPHEUS `ki_n_float` (all n; 78µs; 1e-13). Atkinson F_N moment floor uses E_1, INDEPENDENT of Bickley.
- [phase4_cylinder_peierls.md](phase4_cylinder_peierls.md) — Ki1 vs Ki3, 1/π prefactor, chord branches, rank-N_y white BC, benchmark radii.

### Continuous-µ multi-bounce / specular kernels (Peierls Green's function)
- [phase5_continuous_mu_specular.md](phase5_continuous_mu_specular.md) — Hébert §3.8.3 + Sanchez 2002 Eq.15 = the only textbook continuous-µ multi-bounce forms; Sanchez 1986 TTSP / PS-1982 / Milgram 1978 are the analytical refs.
- [phase5_sanchez_1986_sphere_specular.md](phase5_sanchez_1986_sphere_specular.md) — PDF read. Eq. (A6) IS the continuous-µ multi-bounce specular kernel; T(µ)=1/(1−α·exp(−2aµ)); multi-region extends via piecewise τ(µ).
- [peierls_greens_function_lit.md](peierls_greens_function_lit.md) — Sanchez 2002 §III is LATTICES not sphere; Sanchez 1986 supersedes PS-1982; recommend Variant α (closure-free integral op on Eq. A6); no literature precedent for sphere-specular k_eff via integral-op power iteration.
- [peierls_greens_extensions_lit.md](peierls_greens_extensions_lit.md) — reference list to expand across geom/topology/BC/MG/multi-region; 22 PDFs LOCAL in scratch/literature/; top wishlist Garcia 2021 JCP + Knyazev 1993 + Hébert Ch.6-9.
- [ps1982_and_garcia_extraction.md](ps1982_and_garcia_extraction.md) — PS-1982 Eq.(21) = structurally-indep vacuum-sphere kernel (NO k_eff — use homog-eigval trick + Sood c_crit); Garcia 2021 is fixed-source ONLY (criticality "future work"), gives flux-shape Tables 5/12/18.
- [sanchez_chandrasekhar_three_meanings.md](sanchez_chandrasekhar_three_meanings.md) — disambiguates "transport Green's function": M1 Chandrasekhar angular G / M2 Sanchez source-reduction kernel / M3 trajectory MoC; cylinder M2 needs Carlvik 1965 / Knyazev 1993.

### F_N method / Case eigenfunctions (criticality benchmarks)
- [sood_fn_method_full_extraction.md](sood_fn_method_full_extraction.md) — LA-13511 75-problem catalogue, naming convention, 5-case ramp; F_N machinery lives in KLL 1974 / Westfall-Metcalf / Siewert-Thomas, NOT in the report.
- [sood_2003_vs_1999_extraction.md](sood_2003_vs_1999_extraction.md) — Sood 2003 journal port adds NO method content beyond a general-multigroup k_inf appendix; both versions are TEST SET, not method papers (see L4).
- [sood_2003_cylinder_benchmarks.md](sood_2003_cylinder_benchmarks.md) — `Ua-1-O-CY` (c=1.30) r_c=1.72500292 mfp; `PUb-1-O-CY` (c=1.40) r_c=1.396979 mfp + Φ ratios; Westfall F_N, vacuum-only (no published α∈(0,1) cyl benchmark).
- [kaper_lindeman_leaf_1974_fn_method.md](kaper_lindeman_leaf_1974_fn_method.md) — KLL NSE 54-94 paywalled, NO OA mirror; derivation also in Garcia-Siewert 1979 (paywalled); pivot to Case singular-eigenfunction (B&G Ch.2.6).
- [sphere_fn_method_extraction.md](sphere_fn_method_extraction.md) — Siewert-Thomas 1986: sphere = slab with anti-symmetric BC (one sign flip); F_α recursion + Chebyshev collocation + X-function reusable from slab.
- [meaning1_already_in_fn_method.md](meaning1_already_in_fn_method.md) — KLL X-function + A(ν)/B(ν) Fredholm + `*_angular_flux_from_scalar` ALREADY present for slab+sphere isotropic; gap is API exposure. Check `fn_method/{slab,sphere}/flux_reconstruction.py` first.
- [atkinson_product_nystrom_log_kernels.md](atkinson_product_nystrom_log_kernels.md) — Atkinson 1997 book §4.2 + 1972 Numer.Math. (both LOCAL); Bickley→log via Ki1(τ)~τlogτ; product-Simpson + q=4 graded mesh = the F_N flux-reconstruction fix.

### Anisotropic / reflected criticality (Case-EF & Galerkin reference sets)
- [atalay_1997_reflected_anisotropic.md](atalay_1997_reflected_anisotropic.md) — Atalay 1997 PNE 31(3):229-252, NOT in Zotero. Case singular-EF + Fredholm; UNIFIED slab/sphere via even/odd mode; Tables 2-10 = ready benchmark suite for reflected+anisotropic Sood gap. NOT F_N (see L5) → new `case_method/` package.
- [burkart_ishiguro_siewert_1976_two_region_anisotropic.md](burkart_ishiguro_siewert_1976_two_region_anisotropic.md) — BIS 1976 NSE 61:72-77, PDF read. Case-EF + Chandrasekhar H-function; Tables I (Milne endpoints) + II (critical half-thicknesses), 40 rows each; linear-anisotropy only. NOT F_N.
- [dahl_sjostrand_1979_anisotropic_slab_sphere.md](dahl_sjostrand_1979_anisotropic_slab_sphere.md) — D-S 1979 NSE 69:114-125. Legendre-Galerkin on Carlvik integral eq; full c-eigenvalue SPECTRUM (≤11 eigenvalues) × {µ̄, d}; PRIMARY for P_1 slab+sphere; structurally INDEPENDENT of Sood F_N. NOT F_N.
- [neshat_maiorino_1980_reflected_slab_fn.md](neshat_maiorino_1980_reflected_slab_fn.md) — N-M 1980 ANE 7:79-81, PDF read. Reflected-slab F_N (3 boundary integrals → 3(N+1) system); TRIVIAL extension of bare-slab F_N; unlocks ~10-14 Sood reflected-slab cases; Table 1 critical τ_c values.

### Corpus bibliography (#231 Phase G2)
- [refs_bib_g2_corrections.md](refs_bib_g2_corrections.md) — `docs/refs.bib` COMPLETE (59 entries; 59 rst + 2 docstring-only keys − 2 ruled aliases; pybtex-clean, DOIs CrossRef-verified); ledger of in-text-def errors corrected in FIELDS (keys immutable): WuXieFischer authors, Sanchez2002 coauthors, TTSP vol 15, Stepanek 53-65, Knyazev 368-374, Garcia vols 405/424, BIS 72-77, Carlvik-1965, Mitsis G.J., McCormick-Kuscer + ENDF102 conflations; +MetcalfZweifel1968 = NSE 33(3):318-326 Part II (A19240); swap MUST map PS1982→PomraningSiewert1982, Sood1999→SoodLA13511_1999. Trust the BIB over page/docstring defs.

### Tooling
- [reference_zotero_flakiness.md](reference_zotero_flakiness.md) — diagnose a broken Zotero MCP server (0-hit + conn-refused on port 23119) and fail over to Tier 2 (see L7).

## 3. Active / in-flight state

**No open literature pull.** Most index content is timeless extraction,
not campaign state. The only open SN architecture issue with literature
relevance is **#236** (spatial-scheme ⊗ angular-scheme product) — its
reference basis is [[space_angle_discretization_separability]]. All
other campaign mentions in the topic files (#158/#240, #257/ERR-063,
#196) are CLOSED extraction records, kept for their equation content.

> "Zotero DOWN this session" notes inside topic files are EPHEMERAL
> session state, not durable facts — when an item's annotations matter,
> re-query Zotero; if it answers, capture the annotations.
