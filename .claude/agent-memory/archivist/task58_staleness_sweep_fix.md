# Task #58 — 65-scanner docs-staleness sweep fix (branch `docs/staleness-sweep`)

In-tree edits only (no commit; no pytest; `.claude/`/`scratch/` untouched except
this memo). Sibling of #55/#57 — same per-site behavioral discipline, NO
mechanical blanket rename. 208 findings + 32 false_positives (Haiku scanners,
precision-over-recall); I am the judgment layer — re-verified every finding vs
LIVE tree, rejected ones whose evidence didn't hold.

## GATE OUTPUTS

1. **`.venv/bin/python -m sphinx -W -q docs docs/_build/html`** → **EXIT=0**;
   WARNING/ERROR/CRITICAL count **0** (unchanged from the `-E` pre-edit baseline
   of 0). One transient real WARNING mid-pass (`sood-case-to-geometry` :ref: to a
   title-less anchor — anchors placed AFTER a section title don't resolve as a
   `:ref:` title target; moved the `.. _label:` ABOVE the title → resolved).
2. **`.venv/bin/python -m tests._harness.audit`** → **EXIT=0** (69/69 ERR
   coverage intact; no verifies-linkage touched — matrix.rst NOT regenerated,
   consistent with a pure cross-ref/prose pass).
3. **Mechanical census** (`role_census.py`): **before 120 → after 15 dead
   roles**. All 15 residuals are the verified false-positive classes:
   - **dataclass instance fields / constructor params without defaults** (resolver
     can't getattr-chain them; every one verified live): `SNMethodSpace.quadrature`,
     `ReducedStreamingOperator.requires_upstream_angular_state`, `SNMesh.scheme`,
     `SNMesh.pole_angular_closure`, `SNSolver.inner_solver`, `Solution.angular_flux`,
     `AngularTraceSpace.omega_dot_n`, `HomogeneousResult.representative_energy`,
     `KEigenvalue.eigenvalue_method`, `Mesh1D.widths`.
   - **historical framing** ("no longer exist"): `SNMesh._setup_spherical/_cylindrical`.
   - **planned/future targets**: `orpheus.pn` ("(planned)"), `chord_oracle`
     ("scheduled to promote to"), `derive_second_difference` ("Future: lift...").

## PER-FILE FIX COUNTS (29 .rst + 1 .py)

| File | Fixes | Notes |
|---|---|---|
| collision_probability.rst | 24 | CP test-filename drift (`test_cp_*.py`→`tests/cp/test_*.py`); peierls_nystrom symbol/module fixes (`composite_gl_r`/`lagrange_basis_on_panels` in `.geometry`, `optical_depths_pm` retired→`optical_depth_along_ray`, re-export claim false); `compute_surfaces_1d`→`compute_areas_1d` + `mesh.surfaces`→`mesh.areas`. |
| method_of_characteristics.rst | 22 | short-path `~geometry/numerics/data/sn.*`→`~orpheus.*` (ALL normalized, not just flagged); MoC test-filename drift (`test_moc_*.py`→`tests/moc/test_*.py`); `pwr_pin_equivalent`→`StructuredGeometry.wigner_seitz_pin_cell` (×4). |
| discrete_ordinates.rst | 50→9 | test-path retargets (mms/eigenvalue/phase_c_gates/streaming_operator subdir moves); codec cluster (`solution_to_angular_flux*`/`EquationMap`/`build_equation_map_*` literalized as history); retired sweep-strategy supersession banner (`LegacyTauSymmetricInterpolation`/`BaileyFlatFluxRedist`→MorelMontry sole default); quadrature `Implemented in X`→`Quadrature.factory`; `_mm_weighted_..._single_level`→`compute_psi_half_per_level`; landed-placeholder (`test_phase_d_trajectory_resolvent_crosscheck`); renamed psi_half seed tests. |
| sood_registry.rst | ~14 | THESIS-STALE: `GeometrySpec` class deleted (81b083be)→`geometry_kind: str` + `La13511Case.to_geometry()`; `mixture_to_fn_arrays` retired→direct `Mixture.SigT/SigS/SigP/chi`; "Phase B migration" prospective→accomplished; schema table + bridge-test bullets + rationale rewritten. |
| operator_algebra.rst | ~14 | module fixes (`RadialCharacteristicField`/`DiamondDifference`/`InflowTraceSpace`+`OutflowTraceSpace`→selectors on `AngularTraceSpace`); scattering per-ℓ ladder→`R∘Λ∘M` kernel (`kernel_summands`/`_PerLegendreOrderScattering` retired); malformed `KEigenvalue(A,S,F)`; `build_transport_linear_operator` deleted; never-committed diag paths (repoint `validate_composite_adjoint`→committed test; qualify others "not retained"). |
| singular_eigenfunction.rst | 11 | test-fn name renames (`test_v_se_cyl_*`/`test_v_case_*` suffix drift); `x_function_atalay`→`atalay_X_function` (+wrong module); `test_trajectory_resolvent_cylinder_xverif_sood2003`→`test_peierls_greens_function_cylinder_xverif_sood2003`. |
| galerkin_spectral.rst | 8 | module rename `test_galerkin_spectral_symbolic`→`test_carlvik_galerkin_symbolic` (×8, fn names unchanged). |
| index_convention.rst | 8 | `compute_group_production_rate`/`compute_keff` `:func:`→`:meth:SNSolver.*`; `SNResult`/`SNFixedSourceResult`→`Solution`/`IterationHistory`; `SNSolver.solve_eigenvalue`→`power_iteration` loop; `BulkField.zeros_on`→leaf classes; PR-INDEX-3/PR-INDEX-7 reconciliation (addendum 3). |
| trajectory_resolvent.rst | 6 | `solve_greens_function_specular_sphere`→`_sphere`; 3 `test_trajectory_resolvent_*`→`test_peierls_greens_function_*`; `GeometrySpec` shared-base claim (Billiard/MomentSpace inherit only `object`)→independent math-hearts. |
| fn_method.rst | 5 | test retargets; malformed `sympy.integrate(...)` role; `GeometrySpec`→geometry-tag dispatch; **corrected FALSE claim** "verified against SymPy" → the live test uses `scipy.integrate.quad` singular quadrature. |
| peierls_nystrom.rst | ~8 | `solve_peierls_{cylinder,sphere}_1g`→`geometry.solve_peierls_1g` (CYLINDER_1D/SPHERE_1D dispatch); archived plan/memo paths; `peierls_*.py`→`{slab,cylinder,sphere,geometry}.py`; `_compute_K_bc_sanchez` retired (past-tense the "intentional residual" bullet). |
| thermal_hydraulics.rst | 10 | numbered-dir relics + bare `thermal_hydraulics{,_dae}.py`→`orpheus/thermal_hydraulics/solver{,_dae}.py`. |
| cross_section_data.rst | 4 | `data/`→`orpheus/data/` paths; `gendf.py:281`→`:301`. |
| verification.rst | 2 | eager-tier module subdir paths; MOC `from_annular geometry` (untested claim)→real `test_flux_per_material_matches_scalar`/`test_heterogeneous_flux_depression`. |
| monte_carlo.rst | 3 | `pwr_pin_equivalent`→`wigner_seitz_pin_cell` (×2); `solve_monte_carlo` line range 269-328→589-623 + `monte_carlo.py`→`orpheus/mc/solver.py`. |
| discrete_measures.rst | 7 | `AngularQuadrature` adapters→`Quadrature` class + octants cached-property; section 666-747 "bridge pattern" thesis-rewrite→consolidated `Quadrature` factories + superseded-note. |
| boundary_conditions.rst | 8 | quadrature `GaussLegendre1D`→`Quadrature.gauss_legendre`; malformed `:class:X(args)` roles; `_resolve_bcs`→`realize_boundary_law`; `test_phase_c_gates` subdir; `:func:` module→`:mod:`; historical `AngularQuadrature`/`GaussLegendre1D` literalized. |
| structured_geometry.rst | 3 | `LevelSymmetricSN`/`ProductQuadrature`→`Quadrature.{level_symmetric,product}`; `solve_homogeneous_infinite` module (`derivations.common.eigenvalue`→`homogeneous.solver`). |
| testing/architecture.rst | 2 | `TestZoneSubdivision`/`mesh1d_from_zones`/`TestPWRFactories`/`pwr_pin_equivalent`→live tests + `wigner_seitz_pin_cell`. |
| testing/cross_method.rst | 1 | `mixture_to_fn_arrays`→direct `Mixture` field reads. |
| spherical_harmonics.rst | 2 | test-path subdir; dead `transient-giggling-cake.md` plan link→literalized-as-history (not retained). |
| loss_representations.rst | 2 | `test_loss_action_convention.py` full path; `sn_sweep_strategy.md`→`archive/`. |
| fuel_behaviour.rst | 2 | `data/materials/matpro.py`→`orpheus/`; numbered-dir→`orpheus/fuel/solver.py`. |
| reactor_kinetics.rst | 1 | numbered-dir→`orpheus/kinetics/solver.py`. |
| diffusion_1d.rst | 1 | `kinf_and_spectrum_homogeneous` module (import-site→define-site `common.eigenvalue`). |
| bn_method.rst | 1 | reserved-folder `:mod:` (no `__init__`)→literal (renders dead + prose already "empty"). |
| galerkin_sn_hybrid.rst | 1 | reserved-folder `:mod:`→literal (same). |
| layering.rst | 1 | `AngularFlux.from_flat_with_traces` (moved to `TimedFullField`)→`SNBoundaryRealizer` example + note the codec moved up. |
| peierls.rst | 4 | PDF filename drift (`1982NSE80-481.pdf`, `Sanchez(2002).pdf`); year `Sanchez 1982`→`1977`; `generate_peierls_nystrom_matrix`→`generate_capability_matrices`. |
| **orpheus/sn/operators/streaming.py** | 1 | **#58 core**: D-J retirement note codec roles (`:class:EquationMap`/`:func:build_equation_map_*`/`solution_to_angular_flux*`/`pack_with_traces`)→double-backtick literals (renders via automodule; honest history). |

## REJECTED FINDINGS (evidence did not hold)

| Finding | Why rejected |
|---|---|
| peierls_nystrom.rst:2400 (`peierls_slab`/`peierls_cylinder` `:mod:` "py:mod reference target not found") | The CURRENT clean `-W` build has **no such warning** — bare `:mod:peierls_slab` renders plain-text WITHOUT warning (the finding's warning-evidence is stale/wrong). These are a deliberate legacy-naming convention used ~6× across the page ("reproduces the legacy `peierls_slab`", section titles); 138 is already in the FP list for the same reason. Fixing all 6 = a broad rename beyond flagged scope. LEFT as-is (plain-text, legacy-named, no build impact). |
| sood_registry.rst:687 (`cache` submodule) | false_positive per scanners — verified LIVE (`cache.py` exists, in `__all__`). No action (correctly excluded). |
| discrete_ordinates.rst 6417/11075/15207 (dataclass fields), various dataclass-field/ctor-param census hits | FP class — live fields, resolver getattr-chain limitation. LEFT. |

## AMENDED SUGGESTIONS (overrode a scanner suggestion)

| Site | Scanner said | I wrote instead + why |
|---|---|---|
| discrete_ordinates.rst:5597 / :13450 (addendum + finding) | retarget to `orpheus.numerics.quadrature.AngularQuadrature.spherical_harmonics` | `AngularQuadrature` does NOT exist anywhere (verified `hasattr` False). Correct target is `orpheus.numerics.quadrature.**Quadrature**.spherical_harmonics` (the live class). Scanner conflated the retired Protocol name with the live class. |
| method_of_characteristics.rst:1329 (`pwr_pin_equivalent` "default of 10 fuel + 3 clad + 7 coolant sub-cells") | retarget to `wigner_seitz_pin_cell` | `wigner_seitz_pin_cell` produces a `StructuredGeometry` with region THICKNESSES only — the sub-cell COUNT is a `RegionMesh(n_cells=...)` choice at `Mesh1D.from_geometry`, NOT a method default. Attributing "10/3/7 default" to `wigner_seitz_pin_cell` would be a NEW false claim (caught my own mid-edit error). Wrote: "per-region cell resolution is chosen at `Mesh1D.from_geometry` via `RegionMesh(n_cells=...)`; a typical 10/3/7 subdivision...". |
| fn_method.rst:2007-2008 | retarget test name only | The prose "The closed forms are verified against **SymPy's** `sympy.integrate(...)`" is itself FALSE — the live `test_F_k_primitives_match_scipy_singular_quadrature` verifies against `scipy.integrate.quad` with singularity subdivision. Corrected the CLAIM (SymPy→scipy singular quadrature), not just the role. |
| verification.rst:594 | remove the `from_annular geometry` bullet | Replaced with the two REAL live MOC property tests (`test_flux_per_material_matches_scalar`, `test_heterogeneous_flux_depression`) rather than deleting — the section teaches which properties ARE verified, so a real substitute is better than a gap. |
| operator_algebra.rst 3955-3960 (`kernel_summands`/`_PerLegendreOrderScattering` table row) | remove retired symbols from Code Anchors | The Code-Anchors row (4756/4758) I removed; but the SEPARATE §15.2 table row at 3955 presented the aniso kernel AS `OperatorSum of per-ℓ leaves` (present-tense architecture) — that's thesis-stale, not just a dead anchor. Rewrote to the live `R∘Λ∘M` `OperatorProduct` (Funk-Hecke, retired 93807aa7). |
| operator_algebra.rst 6868 diag path | "commit the diagnostic or cite the equivalent" | The `validate_composite_adjoint` oracle SURVIVED in the committed `tests/sn/operators/test_g_adjoint_reciprocity.py` — repointed there (not "commit the diagnostic", which I can't). |

## SIDE-FLAGS (out of scope; flagged not fixed)

1. **`orpheus/derivations/continuous/sood_registry/la13511.py:311`** — the
   `to_geometry()` infinite-case ERROR MESSAGE string cites
   `orpheus.derivations.common.eigenvalue.solve_homogeneous_infinite`, but that
   function lives at `orpheus.homogeneous.solver`. A `.py` production-code
   staleness (error-message string) — out of my docs scope (only the #58-core
   streaming.py docstring was sanctioned). Numerics/method-implementer fix.
2. **`peierls_slab`/`peierls_cylinder` bare `:mod:` refs** (~6× in
   peierls_nystrom.rst: 138, 716, 1055, 1655, 1928, 2400) — legacy module names
   (live is `peierls_nystrom.{slab,cylinder}`), render plain-text, no build
   warning. A page-wide legacy-naming convention. A dedicated rename pass could
   normalize them, but it's beyond the flagged findings and mostly history-context.
3. **discrete_ordinates.rst leggauss/`build_volume_kernel` in peierls_nystrom.rst §22.9**
   — after `_compute_K_bc_sanchez` deletion, the "residual leggauss consumers"
   list count softened ("The residual call sites..." not "Three call sites"); the
   exact live leggauss-call count in `geometry.py` is uncertain (grep found 0
   direct `leggauss` — possibly aliased import). A follow-up could re-audit the
   exact count.

## BELT-AND-SUSPENDERS grep (found 5 dead roles the scanners MISSED)

A final tree-wide grep for `:func:/:class:/:meth:` roles to the retired symbols I'd
already fixed elsewhere (the census only counts `orpheus.`/`tests.`-prefixed roles,
so short-path `~geometry.*` and thesis-stale present-tense claims slip through)
caught 5 additional sites in UN-FLAGGED sections — all fixed:
- collision_probability.rst:1426 `~geometry.factories.pwr_pin_equivalent` → `wigner_seitz_pin_cell`
- collision_probability.rst:1538 `~geometry.factories.mesh1d_from_zones` → `Mesh1D.from_geometry`
- boundary_conditions.rst:2480, 3134, 920, 941 `SNMesh._resolve_bcs` → `realize_boundary_law`
  (920/941 were present-tense "walks the four faces / builds one AngularTraceSpace" —
  repointed to the shared `resolve_boundary_conditions` driving `realize_boundary_law` per face)
- **singular_eigenfunction.rst:634-648 + code example** — a whole ASPIRATIONAL-design
  section describing a unifying `TransportSolver` Protocol / `GeometrySpec` /
  `from_problem` factory / `geometry_spec` surface as "observed, not posited" — but
  `TransportSolver` and `GeometrySpec` were DELETED in the architectural reset and
  `MomentSpace` has neither `from_problem` nor `geometry_spec` (verified `hasattr` False).
  A code example even `import`ed the deleted `GeometrySpec` module (an UNUSED dead import —
  the example actually uses live `StructuredGeometry`). Fixed: removed the dead import,
  rewrote the "All three classes" block to the honest current shape (`StructuredGeometry` +
  `materials`-dict construction, shared `CriticalSolution`/`FluxSolution` result types — NO
  formal Protocol), + a `.. note::` preserving the "≥3 instances" design rationale as
  history. This was a scope-EXPANSION beyond the flagged findings (scanners caught only the
  trajectory_resolvent sibling + one test name here) — a Cardinal-Rule-1 cluster of false
  present-tense claims. LESSON REINFORCED (→ lessons L-021): the census/scanners miss
  thesis-stale present-tense claims + short-path roles; a targeted tree-wide grep for the
  already-known-retired symbol names is the catch.

## Quality self-assessment (Directive 3)
- Derivation depth 4 · Cross-references 5 (every retarget grep-gated + `hasattr`/
  import-verified against live code; every dead xref repointed or literalized;
  every new `:ref:` anchor confirmed) · Numerical evidence n/a (a cross-ref/
  staleness pass — no flux moves) · Failed approaches 5 (retired-approach history
  preserved via banners+literals, NOT rewritten — quadrature bridge, per-ℓ
  scattering ladder, Phase-B strategy zoo, codec cluster all kept under
  supersession notes) · Code traceability 5 (successors verified live) ·
  Derivation source n/a.
- Weakest: numerical evidence (structural — a staleness pass has no convergence
  table). Same as #55/#57.
