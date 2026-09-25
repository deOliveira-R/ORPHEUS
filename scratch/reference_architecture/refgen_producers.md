# #405 step 2 — the PRODUCER side of reference-solution generation (explorer, 2026-09-24, HEAD 901f64ca)

Probes (all in this scratchpad): `census_types.py` (AST: dataclasses + constructor Call counts),
`reg_dump.py` / `reg_dump2.py` (runtime registry dump: names, phi closure contents, pickle,
materials types), `consumers.py` (AST: registry API call sites), `fieldreads.py` (AST attribute
reads, receiver UNRESOLVED — a homonym ceiling, not a floor), `inline.py` (AST: tests calling
derivations generators directly), `sigs.py` (inspect.signature of entry points).
`[M]` measured; `[R]` inferred.

## 1. Families

| family | entry point(s) | registry | returned type | plain data? |
|---|---|---|---|---|
| homogeneous closed form | `analytical/homogeneous.derive_{1g,2g,4g,2g_n2n}[_continuous]` | `_CASES` (4) AND `_CONTINUOUS` eager (same 4 names) | `VerificationCase` / CRS | VC yes; CRS phi closes over one ndarray |
| flat-source CP (E3/Ki3 CP matrices) | `flat_source_cp/{slab,cylinder,sphere}.all_cases` | `_CASES` (27) | `VerificationCase` | yes (pickle 47/47) |
| SN homogeneous (legacy) | `cases/sn._derive_sn_homogeneous` | `_CASES` (3) | VC | yes |
| MOC / MC homogeneous-cylinder k_inf | `cases/moc`, `cases/mc` | `_CASES` (3 + 9) | VC | yes |
| diffusion sine / transfer matrix | `cases/diffusion.derive_1rg[_continuous]`, `derive_2rg_continuous` | `_CASES` (1) + `_CONTINUOUS` (2) | VC / CRS | CRS phi closes floats/arrays/dicts |
| SN Case singular-eigenfunction (2-region) | `cases/sn.derive_sn_heterogeneous_continuous` | `_CONTINUOUS` (1: `sn_slab_1eg_2rg_S8`) | CRS | closure over arrays/dicts |
| MMS SN (8) + MOC (1) | `mms/sn._build_*_continuous_reference`, `mms/moc._build_moc_mms_continuous_reference`; plus `build_*_mms_case` | `_CONTINUOUS` (9) + inline (102 test calls / 18 files) | CRS wrapping a `*MMSCase` dataclass (13 case types) | NO: closes over the MMSCase, which holds `sigma_t_fn` functions (2 types) and a `Quadrature` |
| Peierls Nyström (unified + native E1) | `peierls_nystrom/cases.continuous_case_builders` (13 thunks), `solve_peierls_{1g,mg}`, `slab.solve_peierls_eigenvalue`, `build_volume_kernel` | lazy `_CONTINUOUS_BUILDERS` (13) + inline (230 test calls / 31 files) | CRS over `PeierlsSolution` / `PeierlsSlabSolution` | PeierlsSolution plain (pickles, 456 B toy) `[M]`; CRS phi a local closure |
| PS-1982 / Atkinson-Nyström | `ps1982_reference`, `fn_method/peierls_atkinson_nystrom` | none | `PS1982Result` | yes |
| trajectory resolvent / Green's fn | `solve_greens_function_{slab,slab_asymmetric,sphere,cylinder,cylinder_mr,annulus,hollow_sphere,…}`, `Billiard` | none; inline (220 test calls / 29 files) | 14 `*Greens*Result` types, `CriticalSolution`/`FluxSolution` via Billiard | yes (frozen arrays) |
| F_N | `solve_fn_slab_bare_critical`, `…sphere…`, `…reflected…`, `compute_kinf_*`, `MomentSpace` | none (capability_rows only); inline 116 calls / 19 files | `SlabFNResult` etc., `CriticalSolution` | yes |
| singular eigenfunction (Case) | `solve_case_method_{slab,sphere}_critical`, `…cylinder_bare_critical`, `Spectrum` | none; inline 49 / 8 | `CaseMethod*Result` | yes |
| Galerkin spectral (Carlvik) | `solve_galerkin_spectral_{slab,sphere}`, `BasisSpace` | none; inline 26 / 6 | `CarlvikGalerkin*Result` | yes |
| Sood LA-13511 + Atalay 1997 | module constants `LA13511_CASES[...]`, `ATALAY_*` | its own dict/tuples (54 La13511Case constructions) | `La13511Case` + `La13511Truth` + sood `Provenance` | yes (Mixture + floats) — published TRUTH, not a generator |
| cross-method catalogue (in tests/) | `tests/gates/cross_method/cases.py` `ALL_CASES`, `case_by_id` | its own list | `CrossMethodCase` | yes |

Registry sizes `[M]` (reg_dump.py): `_CASES` 47 (built 1.8 s); `_CONTINUOUS` eager 16 (1.2 s) + lazy 13; 5 names in both registries (4 homo + `dif_slab_2eg_1rg`), bit-equal k, different `materials` objects.

Traffic `[M]` (consumers.py, AST, positive control = grep's 22 `continuous_get` in 8 files, matched):
registry calls from tests: `get` 76 + `continuous_get` 22 + `continuous_all_names` 4 over 38 files;
direct generator calls from tests (inline.py): ~759 across 8 families. **The registry mediates about 1 in 9 reference acquisitions**; the heavy generators (Peierls, trajectory resolvent) are almost entirely inline.

## 2. Types — field roles

`ContinuousReferenceSolution` (19 constructions, 9 files, 0 in tests):
- problem: `problem: ProblemSpec` {materials, geometry_type, geometry_params, boundary_conditions, is_eigenvalue, n_groups}
- answer: `k_eff`, `phi`, `psi`
- method-of-reference: `tolerance` (prose), `provenance.precision_digits`; MMS smuggles `geometry_params["mms_case"]` (whole problem incl. SN `quadrature`)
- provenance/meta: `name`, `operator_form`, `provenance` {citation, derivation_notes, sympy_expression, precision_digits}, `equation_labels`, `vv_level`, `description`
Readers `[M]` (fieldreads.py, repo-wide, receiver unresolved so an upper bound): `boundary_conditions` 0, `geometry_type` 0, `vv_level` 0 attribute reads (conftest reads it by `getattr` string on `_CASES` objects only), `as_verification_case` 0, `citation`/`derivation_notes`/`sympy_expression`/`precision_digits` 0. In the 8 CRS consumer files: `.problem` 37, `.k_eff` 19, `.geometry_params` 16, `.materials` 13, `.n_groups` 8, `.phi(` 21, `.phi_cell_average` 2, `.phi_on_mesh` 2, `.provenance` 0, `.tolerance` 0, `.psi` 0.

`VerificationCase` (13 constructions, 9 files): problem {method(!), geometry, n_groups, n_regions, materials, geom_params}; answer {k_inf}; meta {name, latex, description, tolerance, vv_level, equation_labels}. `latex`/`description` read only by `generate_rst.py`.

`PeierlsSolution` (frozen; 1 prod + 2 test constructions): answer {r_nodes, phi_values, k_eff, panel_bounds}; method {n_quad_r, n_quad_angular, precision_digits}; problem fragment {cell_radius, geometry_kind, n_groups}. No XS, no BC.

`CylinderGreensMRResult` (frozen, 1 construction): answer {k_eff, psi_g, phi_g, r_nodes, mu_axial_nodes, phi_az_nodes, region_at_node}; method {iterations, converged}. No problem at all (inputs are raw arrays in the call).

`La13511Case` (54): problem {materials, geometry_kind, scattering_order}; answer `truth: La13511Truth` {k_eff_or_kinf, flux_ratios, flux_ratio_groupwise, angular_flux_at_surface, critical_dimension_mfp, extrapolated_endpoint_mfp}; provenance TWICE {sood_table, primary_reference, notes} AND `provenance: Provenance{paper_id, paper_table, primary_reference, notes}`. The geometry SIZE lives in the truth (critical dimension), so `to_geometry()` reads the answer to pose the problem.

`CriticalSolution` / `FluxSolution` (9 / 4 constructions, 0 in tests): the math-heart output; `metadata` dict open-ended.

## Duplications `[M]`

1. **Two `Provenance` classes** with disjoint fields: `common/continuous_reference.py` and `sood_registry/la13511.py` (the sood one re-exported from the package `__init__`).
2. **Geometry kind spelled 7 ways**: ProblemSpec.geometry_type (`sphere-1d`), VerificationCase.geometry (`sph1D`), La13511Case.geometry_kind (`sphere`), PeierlsSolution.geometry_kind, ShippedReference.shape (`sphere-1d`), StructuredGeometry.geometry (`SPH`), CrossMethodCase.geometry (+ `structured_geometry`). Plus bridge maps `_geometry_to_legacy` / `_operator_form_to_legacy_method`.
3. **Literal violations**: ProblemSpec geometry_type values outside `GeometryType`: `sphere`, `cylinder`, `pin-cell-2d` (3 of 8 distinct spellings; 5 construction sites). BC values outside `BoundaryCondition`: `zero_flux` (6 sites), `white_rank2` (2 sites; the deprecated Peierls closure alias, on hollow F.4 refs whose capability row says "vacuum + F.4 cavity closure").
4. **XS spelled 6 ways**: `Mixture` (production type); legacy diffusion dict {transport, absorption, fission, production, chi, scattering} placed raw into `ProblemSpec.materials` (dif_* refs; tests convert with `mixture_from_diffusion_tables`, 6 sites) and re-typed literally in `tests/gates/diffusion/test_properties.py:30`; xs_library `get_xs` dict {sig_t, sig_c, sig_f, nu, chi, sig_s, sig_s1}; raw `sigma_t/sigma_s/nu_sigma_f/chi` arrays (trajectory resolvent); `sig_t/sig_s/nu_sig_f` arrays (Peierls); scalar `c` (F_N, Case, Galerkin); MMS `sigma_t_fn` callables with `materials={}`.
5. **Mixture builders**: 14 in derivations (xs_library.make_mixture + homogeneous `_make_mixture` twin + 3 in mms/sn + 9 `_mix_*` in sood) and 10 test-local builder defs in 8 files.
6. **`_MAT_IDS = {1:[2],2:[2,0],4:[2,3,1,0]}` spelled 6 times** (cases/moc, cases/mc, flat_source_cp ×3, peierls slab); `_RADII` 4 times. `xs_library.get_materials` (0 callers) encodes a DIFFERENT id convention (`len-1-i`: A→3 for 4 regions vs A→2).
7. **Registry enumeration 3 ways**: `_build_registry` hand list (8 modules), `_build_continuous_registry` pkgutil walk (all modules), `tools/verification/generate_capability_matrices._discover_packages` (walks `continuous/*/cases.capability_rows`; 2 packages).

## Selection / metadata `[M]`

- By name only: `get(name)`, `continuous_get(name)`, `LA13511_CASES[case_id]`, `case_by_id`. Name grammars differ per registry (`cp_slab_2eg_2rg`, `sn_mms_*`, `peierls_*` via `naming.reference_name`, Sood `PUa-1-0-SL`).
- Filters: `by_geometry` 0 callers, `by_groups` 0, `by_method` 1 (generate_rst), `continuous_by_operator_form` 0 (and it materialises every lazy Peierls build). No filter on BC, level, tags, cost.
- `vv_level`/`equation_labels` inheritance: `tests/conftest.py` `_resolve_case` resolves `case_name` only through legacy `get`; a test parametrised over a `_CONTINUOUS`-only name inherits nothing.
- Dead/stale: `continuous_register` 0 callers, docstring names a retired `_continuous_modules` list (reference_values.py:296); `mms/sn.all_cases` (line 241) 0 consumers and not walked by `_build_registry`; `as_verification_case` 0 callers.

## Q-R1 relevance

- Four math-heart classes already take the method-independent posing `(StructuredGeometry, dict[int, Mixture], method knobs)`: `MomentSpace`, `Billiard`, `Spectrum`, `BasisSpace` (sigs.py). That is the Q-R1 shape one layer up from `MaterialMesh`.
- `MaterialMesh` appears in 0 files under `orpheus/derivations` (grep; positive control: 15 files elsewhere). No reference emits a mesh; `sood_registry.builders.build_mesh` emits `Mesh1D`.
- MMS references own the METHOD head: 13 of 13 MMS case types carry discretisation choices (`quadrature` in the 12 SN types; `n_azi/n_polar/ray_spacing` in MOC), contrary to the Q-R1 split.
- `_CASES` 47/47 pickle; `_CONTINUOUS` eager 0/16 pickle (local-closure phi).
