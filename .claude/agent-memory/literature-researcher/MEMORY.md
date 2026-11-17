# Literature Researcher — Memory Index

Section 2 is a negative-lookup inventory: before extracting a paper, search this
page for it. Each row links one topic file; the link text names the subject and
the hook says what the file settles. The content is in the file, never here.
A paper filed under another paper's name is listed in section 3.

## 1. Lessons

- [lessons.md](lessons.md) — Eight lessons no carrier states, and a ledger of where the other parts are carried.

## 2. Extraction inventory by topic

### Curvilinear and spherical SN
- [The "Bailey 2009, Ann. Nucl. Energy 35" citation is refuted](bailey_ane_chimera_citation_refuted.md) — A fusion of three records; its equations are Bailey-Morel-Chang 2010 Eqs. 50, 52, 74.
- [Lathrop 2000, the angular-scheme taxonomy](lathrop_2000_angular_scheme_taxonomy.md) — Five sphere schemes; any τ≠½ is first-order in angle; also Chaland-Samba 2016, Larsen-Morel §1.5.1.
- [Morel-Montry τ and angular cell edges](morel_montry_tau_angular_cell_edges.md) — Sphere and cylinder τ rulings from Morel-Montry 1984, Bailey-Morel-Chang 2010 and Hébert.
- [Reed-Lathrop 1970, angular truncation error](reed_lathrop_1970_angular_truncation.md) — Their weighted diamond solves for the ordinates; cylinder form; two solve-free gates.
- [Linear discontinuous, slab and curvilinear](issue_158_linear_discontinuous_cell_update.md) — Larsen-Morel 1989 slab; Adams-Martin 1992 App. A spherical; Palmer-Adams 1993 diffusion limit.
- [Multi-dimensional linear discontinuous closure](multi_d_ld_closure.md) — Bilinear or trilinear discontinuous elements, not simplex P1, which fails the thick-diffusion limit.
- [Sphere SN pole closure](sphere_sn_pole_closure_canonical.md) — Hébert §3.9.4, Eqs. 3.418–3.439, is the canonical spherical closure.
- [Carlson coupled pole](phase_d_carlson_coupled_pole.md) — Hébert Eqs. 3.432–3.435 verbatim; the μ=−1 seed; cylinder analogue 3.407–3.410.
- [Curvilinear sweep directness](curvilinear_sweep_directness_ruling.md) — Ruling: the sweep is direct, one pass; the only published lag is the albedo boundary.
- [Central cell at r=0](spherical_sn_central_cell_spatial_order.md) — Hébert and Stacey both use plain midpoint diamond; Eq. 3.430 flux is a volume average.
- [Spatial order at the origin](sphere_sn_spatial_order_at_origin.md) — No canonical second-order central-cell spatial closure; the lead is Wu-Xie-Fischer 1999.
- [Space and angle separability](space_angle_discretization_separability.md) — Spatial (Larsen-Morel-Miller 1987) and angular (Bailey-Morel-Chang 2010) diffusion-limit conditions are separate axes.
- [Face transmission as a Padé approximant](face_transmission_pade.md) — Lesaint-Raviart 1974 Theorem 2; Morel-Wareing-Smith 1996 Eq. 74; diamond difference is Padé [1/1].

### Angular quadrature
- [Galerkin quadrature and the linear-discontinuous mass matrix](galerkin_quadrature_metric_and_ld_mass.md) — Morel 1989 and 2017, Shands 2025, Sanchez-Ragusa 2011; Adams-Martin (A.4) mass integrals.
- [Tabuchi-Yamamoto polar quadrature](ty_polar_quadrature_moc.md) — Yamamoto 2007 JNST 44(2) is not local; table provenance, weight convention, access status.
- [LA-3186 and LA-4058, level-symmetric quadrature](la3186_level_symmetric_quadrature.md) — Axis-weight construction, Table I contents, point-weight positivity up to S20, printed rounding slips.

### SN solver and sweep theory
- [Diffusion synthetic acceleration, seven primaries](dsa_primaries_full_extraction.md) — Alcouffe 1977, Larsen 1982, McCoy-Larsen 1982, Morel 1982, Adams-Martin 1992, Warsa 2004; memo path.
- [Adams and Larsen 2002](adams_larsen_2002_iteration_theory.md) — Where each convergence result sits in the review; notation traps; memo path.
- [Morel 1989, sweep versus linear system](morel_1989_si_vs_apply_equivalence.md) — Page 75: source iteration equals the linear system only for lower-triangular angular leakage.
- [Larsen and Morel 2010 review](larsen_morel_2010_review_extraction.md) — Equation map; no cylinder SN, adjoint or verification content; printed-page offset.
- [Method naming evidence](method_canonical_naming_evidence.md) — Verbatim quotes for the "singular eigenfunction expansion" and "trajectory-based" method names.

### Nuclear data formats
- [ENDF-6, GENDF and NJOY: the (n,2n) rulings](endf6_gendf_njoy_n2n_formats.md) — Manual URLs and page offsets; Legendre convention; GENDF record layout; nine-nuclide ENDF/B-VIII.0 census.

### Time (α) eigenvalue
- [α eigenvalue existence and spectrum](alpha_eigenvalue_existence_spectrum.md) — Lehner-Wing: the slab α₀ always exists; where non-existence arises; Dorning 2010 mis-states it.

### Multigroup sources and collapse
- [Fission-source χ indexing](peierls_mg_fission_source_chi_indexing.md) — Hébert Eqs. 3.57–3.58: the emission spectrum is local to the source point.
- [Condensation formulas](energy_condensation_collapse_formulas.md) — Hébert §3.5 and Stamm'ler Ch. VI formulas; group nesting; NJOY weighting; Rahnema 2008 projection.
- [Adjoint bilinear collapse](adjoint_bilinear_collapse_p6.md) — Bell-Glasstone Ch. 6, Eqs. 6.135–6.136 bilinear weighting; fission dyad; perturbation grounding.
- [Spherical-harmonic moment projection](sh_flux_moment_projection_root_cause.md) — The moment matrix exists because scattering is anisotropic; Hébert, Brockmann, Ahrens, Fletcher equations.

### Interface currents
- [Sanchez-McCormick 1982, rank-N currents](rank_n_interface_current_canonical.md) — §III.F is the only rank-N interface-current source in the textbook corpus.
- [Sanchez-McCormick 1982, Eq. 166](sanchez_mccormick_rank_n_per_face.md) — The Eq. 166 cluster verbatim; the N=2 residual is a missing μ weight.
- [Rank-N closure in four references](rank_n_closure_four_references_synthesis.md) — Ligou 1982, Sanchez 2002, Stamm'ler Ch. IV and Stacey 2007 are all scalar.
- [Hébert 2009 Ch. 3, interface currents](hebert_2009_ch3_interface_currents.md) — Angular expansion only for 2-D Cartesian; 1-D curvilinear is scalar; Villarino-Stamm'ler renormalisation.
- [Stamm'ler 1983, Ch. VI and Ch. IV](stammler_1983_ch6_interface_currents.md) — Ch. VI is SN; Ch. IV is collision probability with scalar currents.
- [Sanchez 1977, NSE 64](sanchez_1977_nse64_canonical_ic.md) — Rank-3 surface modes, fixed source only; accuracy limits in Tables VI–VII.
- [No curvilinear double-P_N reference tables](dpn_curvilinear_nonexistence.md) — None exist for cylinder or sphere; Sood 2003 is the surrogate.
- [Rank-N curvilinear leads](rank_n_ic_curvilinear_literature_leads.md) — Bogado Leite 1998, Krishnani 1982, Mohanakrishnan 1982: read before any novelty claim.
- [Lambert and Marshak bases](direction_q_lambert_marshak_derivation.md) — Sanchez 2014 NSE 177 is differential P_N, not integral transport; the change-of-basis derivation.
- [Knyazev 1993](knyazev_1993_cylinder_anisotropic_ic.md) — Atomic Energy 74(5):368–374, the cylinder Ki_{2+k} source; "Knyazev-Selivanov 2014" does not exist.

### Collision probability
- [Collision-probability moment integrals](cp_moment_integrals.md) — Slab and cylinder closed forms from Abramowitz-Stegun §5.1 and Hébert §3.4–3.5.
- [Carlvik 1967](carlvik_1967_finite_cylinder_cp.md) — A two-page note with E_3 kernels; the cylindrical-cell paper is Carlvik 1965, Geneva.
- [Sphere collision probability and the white boundary](bickley_naylor_sphere_white_bc.md) — The sphere kernel is E_2, not a Bickley function; Hébert §3.8.5.
- [Bickley function libraries](bickley_function_libraries.md) — Davierwalla 1982, Lorensi 2025 (Ki₁ and Ki₂ only) and Amos 1983 against the in-tree functions.
- [Cylinder Peierls equation](phase4_cylinder_peierls.md) — Ki₁ against Ki₃, the 1/π prefactor, the two chord branches, the white boundary.

### Continuous-μ specular reflection
- [Continuous-μ specular kernels](phase5_continuous_mu_specular.md) — Hébert §3.8.3 and Sanchez 2002 Eq. 15 multi-bounce forms; the analytical references.
- [Sanchez 1986, sphere](phase5_sanchez_1986_sphere_specular.md) — Eq. A6 is the continuous-μ multi-bounce specular kernel; the multi-region extension.
- [Green's-function literature for the sphere](peierls_greens_function_lit.md) — Sanchez 2002 §III treats lattices, not spheres; the recommended formulation rests on Sanchez 1986.
- [Extensions wishlist](peierls_greens_extensions_lit.md) — Candidate extensions by geometry and boundary condition; which references are local, which wanted.
- [Pomraning-Siewert 1982 and Garcia 2021](ps1982_and_garcia_extraction.md) — Pomraning-Siewert Eq. 21 vacuum-sphere kernel; Garcia 2021 fixed-source Tables 5, 12, 18.
- [Three meanings of "Green's function"](sanchez_chandrasekhar_three_meanings.md) — Chandrasekhar's angular function, Sanchez's source-reduction kernel, the trajectory method.

### F_N method and criticality benchmarks
- [Sood LA-13511 (1999)](sood_fn_method_full_extraction.md) — The 75-problem catalogue and its naming convention; the method is in its references.
- [Sood 2003 against 1999](sood_2003_vs_1999_extraction.md) — The journal version adds only a multigroup k_inf appendix.
- [Sood 2003, cylinder values](sood_2003_cylinder_benchmarks.md) — Critical radii for the two one-group cylinder problems; vacuum boundary only.
- [Kaper-Lindeman-Leaf 1974](kaper_lindeman_leaf_1974_fn_method.md) — Paywalled with no open copy; where else the F_N derivation is printed.
- [Siewert-Thomas 1986, sphere](sphere_fn_method_extraction.md) — The sphere is the slab with an antisymmetric boundary condition; collocation details.
- [F_N angular-flux reconstruction](meaning1_already_in_fn_method.md) — The reconstruction already exists in the tree's `fn_method` package; the gap is exposure.
- [Atkinson product-Nyström](atkinson_product_nystrom_log_kernels.md) — Atkinson 1997 §4.2 and 1972: product integration on a graded mesh for logarithmic kernels.

### Anisotropic and reflected criticality
- [Atalay 1997](atalay_1997_reflected_anisotropic.md) — Prog. Nucl. Energy 31(3):229–252; Case singular eigenfunctions, not F_N; Tables 2–10.
- [Burkart-Ishiguro-Siewert 1976](burkart_ishiguro_siewert_1976_two_region_anisotropic.md) — NSE 61:72–77; Case eigenfunctions, not F_N; Tables I and II, 40 rows each.
- [Dahl-Sjöstrand 1979](dahl_sjostrand_1979_anisotropic_slab_sphere.md) — NSE 69:114–125; Legendre-Galerkin eigenvalue spectrum for the linearly anisotropic slab and sphere.
- [Neshat-Maiorino 1980](neshat_maiorino_1980_reflected_slab_fn.md) — Ann. Nucl. Energy 7:79–81; reflected-slab F_N; Table 1 critical thicknesses.

### Bibliography and mathematical foundations
- [Corrections ledger for `docs/refs.bib`](refs_bib_g2_corrections.md) — Per-entry corrections; trust the bib over an in-text definition; two key consolidations.
- [Riesz, adjoint, pseudoinverse and frame-letter citations](math_foundations_citations_riesz_penrose.md) — Kreyszig, Lee, Penrose and Moore page citations; which authors write T for synthesis.
- [Frame and kernel hierarchy sources](frame_kernel_hierarchy_sources.md) — Frame, Markov-kernel and coupling definitions read across two passes; page offsets.

### Verification sources
- [Numerical-error taxonomy sources](verification_error_ontology_sources.md) — Oberkampf-Trucano 2002, Arioli 2013, Atkinson 1997, Trefethen 2008 locators; the paywalled list.
- [Verification-case anatomy sources](verification_case_anatomy_sources.md) — Oberkampf-Trucano 2007, ICSBEP, Ganapol 2008, ANL-7416 identifiers, build-system trace papers.

### Tooling
- [The full NSE run is on disk as zip bundles](user_nse_volume_archive.md) — Search the zip listings before calling an NSE paper not local; file names are truncated.
- [Zotero when the server is down](reference_zotero_flakiness.md) — A read-only fallback through a scratchpad copy of `zotero.sqlite`.

## 3. Papers filed under another name

- Chaland-Samba 2016 (two-angle closure), Hu-Azmy 2020–2021 and Wu-Xie-Fischer 1999: `lathrop_2000_angular_scheme_taxonomy.md`.
- Adams-Martin 1992 App. A, Palmer-Adams 1993 and Morel-Wareing-Smith 1996: `issue_158_linear_discontinuous_cell_update.md`.
- The Bailey-Morel-Chang 2010 equation map and the cylinder cosine map: `bailey_ane_chimera_citation_refuted.md`.
- Westfall-Metcalf 1973 and Mitsis 1963, as naming quotes: `method_canonical_naming_evidence.md`.
- Boyd 2019 (OpenMC multiplicity), ISOTXS and the WIMS library format: `endf6_gendf_njoy_n2n_formats.md`.
- Christensen 2016, Casazza 2000 and Casazza-Lynch 2016, the operator letters: `math_foundations_citations_riesz_penrose.md`.

## 4. Primaries named and still unread

- Grant 1968, J. Comput. Phys. 2(4):381–402, origin of the weighted-diamond node: `reed_lathrop_1970_angular_truncation.md`.
- Dudziak-O'Dell-Alcouffe 1979, LA-7911-PR, the starting-direction primary: `morel_montry_tau_angular_cell_edges.md`.
- Bailey-Adams-Chang 2008, LLNL-CONF-407632, a spatial paper: `bailey_ane_chimera_citation_refuted.md`.
- Bogado Leite 1998, rank-N currents in curvilinear geometry: `rank_n_ic_curvilinear_literature_leads.md`.
- The body of Yamamoto 2007 (fit order, range and bound): `ty_polar_quadrature_moc.md`.

## 5. Standing notes

- An issue's state comes from `gh issue view`, and a campaign's from git, never from a topic file.
- A "Zotero down" note in a topic file describes one session; query Zotero again.
