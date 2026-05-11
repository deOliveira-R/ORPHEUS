Verification Matrix
===================

.. note::

   Auto-generated from ``tests._harness.registry.TEST_REGISTRY``
   by ``tools/verification/generate_matrix.py``. Do not edit by
   hand — changes will be overwritten on the next rebuild.

Total tests collected: **3402**

V&V level distribution
----------------------

.. csv-table::
   :header: Level, Count, Share
   :widths: 15, 10, 10

   L0, 900, 26.5%
   L1, 837, 24.6%
   L2, 36, 1.1%
   L3, 0, 0.0%
   foundation, 1601, 47.1%
   unmarked, 28, 0.8%

Tagging source
--------------

How each test acquired its V&V level (see ``tests/conftest.py`` for the precedence chain).

.. csv-table::
   :header: Source, Count
   :widths: 20, 10

   explicit, 3295
   verify, 0
   class-name, 46
   func-name, 0
   case, 33
   unmarked, 28

Module × level grid
-------------------

.. csv-table::
   :header: Module, L0, L1, L2, L3, FD, ??
   :widths: 40, 6, 6, 6, 6, 6, 6

   cp/test_boundary_realizer_stub, 3, 0, 0, 0, 0, 0
   cp/test_cylinder, 0, 9, 0, 0, 0, 0
   cp/test_cylinder_pss, 0, 0, 0, 0, 16, 0
   cp/test_diagnostics, 8, 28, 0, 0, 0, 0
   cp/test_peierls_cylinder_flux, 0, 4, 0, 0, 0, 0
   cp/test_peierls_flux, 0, 1, 0, 0, 0, 0
   cp/test_peierls_rank_n_protocol, 0, 20, 0, 0, 0, 0
   cp/test_peierls_sphere_flux, 0, 4, 0, 0, 0, 0
   cp/test_properties, 12, 0, 0, 0, 0, 0
   cp/test_slab, 0, 9, 0, 0, 0, 0
   cp/test_sphere, 0, 9, 0, 0, 0, 0
   cp/test_verification, 4, 25, 5, 0, 0, 0
   cross_method/test_eigenvalue, 0, 31, 0, 0, 53, 0
   cross_method/test_polymorphism, 0, 0, 0, 0, 5, 0
   data/test_cross_section_data, 11, 0, 0, 0, 0, 0
   data/test_mixture, 4, 0, 0, 0, 0, 0
   data/test_mixture_scattering_ratio, 0, 0, 0, 0, 3, 0
   derivations/test_atkinson_product_nystrom, 0, 5, 0, 0, 3, 3
   derivations/test_capability_matrices, 0, 0, 0, 0, 3, 0
   derivations/test_carlvik_galerkin_slab, 0, 18, 0, 0, 0, 0
   derivations/test_carlvik_galerkin_sood_registry, 0, 5, 0, 0, 2, 0
   derivations/test_carlvik_galerkin_sphere, 0, 18, 0, 0, 0, 0
   derivations/test_carlvik_galerkin_symbolic, 0, 0, 0, 0, 8, 0
   derivations/test_carlvik_galerkin_xverif_fn, 0, 10, 0, 0, 0, 0
   derivations/test_case_method_slab, 0, 17, 0, 0, 1, 0
   derivations/test_case_method_slab_sphere_parity_flip, 0, 4, 0, 0, 0, 0
   derivations/test_case_method_sphere, 0, 4, 0, 0, 1, 0
   derivations/test_case_method_symbolic, 0, 0, 0, 0, 9, 0
   derivations/test_case_method_x_function, 2, 3, 0, 0, 0, 0
   derivations/test_case_method_z0, 0, 11, 0, 0, 0, 0
   derivations/test_cp_geometry, 48, 0, 0, 0, 0, 0
   derivations/test_fn_la13511_kinf, 0, 0, 0, 0, 17, 0
   derivations/test_fn_la13511_slab, 0, 0, 0, 0, 21, 0
   derivations/test_fn_la13511_slab_flux, 0, 10, 0, 0, 0, 0
   derivations/test_fn_la13511_slab_flux_symbolic, 0, 0, 0, 0, 6, 0
   derivations/test_fn_la13511_slab_reflected, 0, 18, 0, 0, 4, 1
   derivations/test_fn_la13511_slab_xverif, 0, 2, 0, 0, 0, 0
   derivations/test_fn_la13511_sphere, 0, 0, 0, 0, 11, 0
   derivations/test_fn_la13511_sphere_flux, 0, 10, 0, 0, 0, 0
   derivations/test_fn_la13511_sphere_xverif, 0, 3, 0, 0, 0, 0
   derivations/test_fn_method_moment_space, 0, 0, 0, 0, 14, 0
   derivations/test_fn_projection_vs_kll_flux, 0, 4, 0, 0, 4, 2
   derivations/test_fn_sood_table10_symmetric_pu_h2o, 0, 1, 0, 0, 2, 0
   derivations/test_galerkin_spectral_basis_space, 0, 0, 0, 0, 18, 0
   derivations/test_kernels, 55, 0, 0, 0, 0, 0
   derivations/test_la13511_to_geometry, 0, 0, 0, 0, 131, 0
   derivations/test_mu_weighted_basis, 0, 0, 0, 0, 1, 0
   derivations/test_path_ai_legacy_plain_gl_signature, 0, 3, 0, 0, 1, 0
   derivations/test_peierls_assembly_drivers, 0, 0, 0, 0, 9, 0
   derivations/test_peierls_closure_operator, 0, 0, 0, 0, 48, 0
   derivations/test_peierls_convergence, 5, 0, 0, 0, 0, 0
   derivations/test_peierls_cylinder_eigenvalue, 3, 5, 0, 0, 0, 0
   derivations/test_peierls_cylinder_g_bc_3d_symbolic, 0, 5, 0, 0, 0, 0
   derivations/test_peierls_cylinder_geometry, 10, 0, 0, 0, 0, 0
   derivations/test_peierls_cylinder_knyazev_symbolic, 0, 22, 0, 0, 0, 0
   derivations/test_peierls_cylinder_multi_region, 7, 0, 0, 0, 3, 0
   derivations/test_peierls_cylinder_prefactor, 4, 0, 0, 0, 0, 0
   derivations/test_peierls_cylinder_white_bc, 4, 3, 0, 0, 4, 0
   derivations/test_peierls_geometry, 0, 0, 0, 0, 32, 0
   derivations/test_peierls_greens_function_annulus_solver, 0, 13, 0, 0, 0, 0
   derivations/test_peierls_greens_function_annulus_symbolic, 0, 0, 0, 0, 22, 0
   derivations/test_peierls_greens_function_cylinder_solver, 0, 13, 0, 0, 0, 0
   derivations/test_peierls_greens_function_cylinder_symbolic, 0, 0, 0, 0, 9, 0
   derivations/test_peierls_greens_function_cylinder_xverif_sood2003, 0, 2, 0, 0, 0, 0
   derivations/test_peierls_greens_function_garcia2021, 0, 0, 0, 0, 17, 0
   derivations/test_peierls_greens_function_hollow_sphere_solver, 0, 14, 0, 0, 0, 0
   derivations/test_peierls_greens_function_hollow_sphere_symbolic, 0, 0, 0, 0, 18, 0
   derivations/test_peierls_greens_function_mg, 0, 0, 0, 0, 7, 0
   derivations/test_peierls_greens_function_mr, 0, 0, 0, 0, 4, 0
   derivations/test_peierls_greens_function_slab_asymmetric_solver, 0, 14, 0, 0, 0, 0
   derivations/test_peierls_greens_function_slab_asymmetric_symbolic, 0, 0, 0, 0, 16, 0
   derivations/test_peierls_greens_function_slab_solver, 0, 14, 0, 0, 0, 0
   derivations/test_peierls_greens_function_slab_symbolic, 0, 0, 0, 0, 10, 0
   derivations/test_peierls_greens_function_solver, 0, 1, 0, 0, 3, 0
   derivations/test_peierls_greens_function_symbolic, 0, 0, 0, 0, 9, 0
   derivations/test_peierls_greens_function_vacuum, 0, 0, 0, 0, 5, 0
   derivations/test_peierls_greens_function_xverif, 0, 5, 0, 0, 3, 0
   derivations/test_peierls_greens_function_xverif_ps1982, 0, 0, 0, 0, 6, 0
   derivations/test_peierls_multigroup, 9, 10, 0, 0, 8, 0
   derivations/test_peierls_nystrom_verification, 0, 4, 0, 0, 0, 0
   derivations/test_peierls_rank2_bc, 9, 24, 0, 0, 4, 0
   derivations/test_peierls_rank_n_bc, 59, 15, 0, 0, 0, 0
   derivations/test_peierls_rank_n_class_b_mr_mg, 0, 31, 0, 0, 0, 0
   derivations/test_peierls_rank_n_conservation, 0, 0, 0, 0, 4, 0
   derivations/test_peierls_rank_n_primitives, 19, 0, 0, 0, 0, 0
   derivations/test_peierls_reference, 40, 25, 0, 0, 16, 0
   derivations/test_peierls_slab_legacy_aggregate, 0, 0, 0, 0, 4, 0
   derivations/test_peierls_specular_bc, 0, 0, 0, 0, 27, 0
   derivations/test_peierls_specular_continuous_mu_symbolic, 0, 0, 0, 0, 4, 0
   derivations/test_peierls_specular_slab_symbolic, 0, 10, 0, 0, 0, 0
   derivations/test_peierls_specular_symbolic, 0, 18, 0, 0, 0, 0
   derivations/test_peierls_sphere_eigenvalue, 0, 4, 0, 0, 0, 0
   derivations/test_peierls_sphere_geometry, 21, 0, 0, 0, 0, 0
   derivations/test_peierls_sphere_prefactor, 6, 0, 0, 0, 0, 0
   derivations/test_peierls_sphere_white_bc, 0, 4, 0, 0, 0, 0
   derivations/test_peierls_variant_alpha_core, 0, 0, 0, 0, 8, 0
   derivations/test_quadrature, 7, 0, 0, 0, 44, 0
   derivations/test_singular_eigenfunction_cylinder, 0, 8, 0, 0, 14, 0
   derivations/test_singular_eigenfunction_cylinder_xverif, 0, 1, 0, 0, 0, 0
   derivations/test_singular_eigenfunction_spectrum, 0, 0, 0, 0, 16, 0
   derivations/test_sn_mms_anisotropic_symbolic, 0, 0, 0, 0, 12, 0
   derivations/test_sood_registry_cache, 0, 0, 0, 0, 15, 0
   derivations/test_sood_registry_compatibility, 0, 2, 0, 0, 107, 0
   derivations/test_sood_registry_wide_bare_critical, 0, 17, 0, 0, 2, 0
   derivations/test_sood_registry_wide_kinf, 0, 0, 0, 0, 49, 0
   derivations/test_trajectory_resolvent_billiard, 0, 0, 0, 0, 11, 0
   derivations/test_trajectory_resolvent_chord_oracle, 0, 0, 0, 0, 21, 0
   derivations/test_trajectory_resolvent_power_iterate, 0, 0, 0, 0, 6, 0
   diffusion/test_boundary_realizer_stub, 3, 0, 0, 0, 0, 0
   diffusion/test_continuous_reference, 0, 8, 0, 0, 0, 0
   diffusion/test_diffusion, 0, 2, 0, 0, 0, 0
   diffusion/test_properties, 3, 0, 0, 0, 0, 0
   geometry/test_bc_equivalence_snapshot, 0, 15, 0, 0, 0, 0
   geometry/test_bc_errors, 0, 0, 0, 0, 11, 0
   geometry/test_bc_universal_invariants, 0, 30, 0, 0, 11, 0
   geometry/test_bound_compat, 10, 0, 0, 0, 0, 0
   geometry/test_boundary, 0, 0, 0, 0, 24, 0
   geometry/test_boundary_trace_law, 0, 0, 0, 0, 16, 0
   geometry/test_geometry, 0, 0, 0, 0, 55, 0
   geometry/test_mesh, 0, 0, 0, 0, 10, 0
   geometry/test_reduced_operator, 0, 0, 0, 0, 49, 0
   geometry/test_structured_geometry, 0, 0, 0, 0, 36, 0
   homogeneous/test_continuous_reference, 0, 7, 0, 0, 0, 0
   homogeneous/test_homogeneous, 0, 4, 0, 0, 0, 0
   l1_analytical/test_kinf_homogeneous, 0, 15, 0, 0, 0, 0
   l1_analytical/test_mms_curvilinear_aniso_dd_convergence, 0, 2, 0, 0, 0, 0
   l1_analytical/test_pole_closure_flat_flux_identity, 0, 5, 0, 0, 0, 0
   mc/test_boundary_realizer_stub, 3, 0, 0, 0, 0, 0
   mc/test_convergence, 0, 0, 3, 0, 0, 0
   mc/test_cross_verification, 0, 0, 2, 0, 0, 0
   mc/test_gaps, 7, 9, 0, 0, 0, 0
   mc/test_monte_carlo, 0, 12, 0, 0, 0, 0
   mc/test_properties, 24, 0, 0, 0, 0, 0
   moc/test_boundary_realizer_stub, 3, 0, 0, 0, 0, 0
   moc/test_mms, 0, 3, 0, 0, 0, 0
   moc/test_moc, 0, 3, 0, 0, 0, 0
   moc/test_properties, 4, 0, 0, 0, 0, 0
   moc/test_quadrature, 24, 0, 0, 0, 0, 0
   moc/test_ray_tracing, 22, 0, 0, 0, 0, 0
   moc/test_verification, 27, 15, 6, 0, 0, 0
   numerics/test_diagonal_operator, 17, 0, 0, 0, 0, 0
   numerics/test_incoming_ordinate_mask_tensor, 13, 0, 0, 0, 0, 0
   numerics/test_iteration, 0, 1, 0, 0, 10, 0
   numerics/test_measure, 0, 16, 0, 0, 32, 0
   numerics/test_measure_partition, 12, 0, 0, 0, 0, 0
   numerics/test_operator, 0, 0, 0, 0, 55, 0
   numerics/test_periodic_wrap_operator, 5, 0, 0, 0, 0, 0
   numerics/test_permutation_operator, 9, 1, 0, 0, 0, 0
   numerics/test_projection_operators, 16, 6, 0, 0, 0, 0
   numerics/test_registry, 0, 0, 0, 0, 37, 0
   numerics/test_registry_mixin, 0, 0, 0, 0, 10, 0
   numerics/test_rules_1d, 0, 5, 0, 0, 16, 0
   numerics/test_rules_product, 0, 3, 0, 0, 14, 0
   numerics/test_rules_sphere, 0, 7, 0, 0, 29, 0
   numerics/test_space, 0, 0, 0, 0, 17, 0
   numerics/test_spherical_harmonics, 5, 7, 0, 0, 0, 0
   numerics/test_symmetry, 0, 0, 0, 0, 71, 0
   numerics/test_tensor_product_operator, 18, 0, 0, 0, 0, 0
   numerics/test_trace_space, 13, 7, 0, 0, 0, 0
   regression/test_dd_regression, 0, 0, 0, 0, 0, 11
   sn/test_2d_octant_sweep_equivalence, 0, 7, 0, 0, 0, 0
   sn/test_angular_average_operator, 13, 4, 0, 0, 0, 0
   sn/test_boundary_conditions, 0, 0, 0, 0, 0, 11
   sn/test_boundary_realize, 0, 3, 0, 0, 8, 0
   sn/test_cartesian, 1, 6, 0, 0, 0, 0
   sn/test_cell_update_batch, 10, 0, 0, 0, 0, 0
   sn/test_cylindrical, 4, 10, 11, 0, 0, 0
   sn/test_discrete_ordinates_2d, 0, 0, 2, 0, 0, 0
   sn/test_fission_operator, 0, 0, 0, 0, 10, 0
   sn/test_heterogeneous_transport, 0, 2, 0, 0, 0, 0
   sn/test_legendre_moment_scattering, 9, 0, 0, 0, 0, 0
   sn/test_method_space, 6, 0, 0, 0, 0, 0
   sn/test_mms, 0, 2, 0, 0, 0, 0
   sn/test_mms_2d, 0, 3, 0, 0, 0, 0
   sn/test_mms_aniso, 0, 2, 0, 0, 0, 0
   sn/test_mms_curvilinear, 0, 2, 0, 0, 0, 0
   sn/test_mms_heterogeneous, 0, 4, 0, 0, 0, 0
   sn/test_octants_property, 60, 0, 0, 0, 0, 0
   sn/test_properties, 4, 0, 0, 0, 0, 0
   sn/test_quadrature, 49, 0, 0, 0, 0, 0
   sn/test_scattering_operator, 0, 0, 0, 0, 17, 0
   sn/test_sn_boundary_realizer, 0, 19, 0, 0, 0, 0
   sn/test_snmesh_consumes_reduced, 0, 0, 0, 0, 15, 0
   sn/test_snmesh_realizer_wiring, 0, 9, 0, 0, 0, 0
   sn/test_snmesh_sweep_graphs, 13, 0, 0, 0, 0, 0
   sn/test_snstreamingoperator, 0, 0, 0, 0, 30, 0
   sn/test_solver_components, 41, 0, 0, 0, 0, 0
   sn/test_spherical, 13, 7, 6, 0, 0, 0
   sn/test_sweep_graph, 63, 0, 0, 0, 0, 0
   sn/test_sweep_operator_inconsistency, 0, 4, 0, 0, 0, 0
   sn/test_sweep_regression, 12, 0, 0, 0, 0, 0
   sn/test_unified_sweep_dispatch, 0, 0, 0, 0, 9, 0
   spatial/test_boundary_face_flux, 0, 0, 0, 0, 21, 0
   spatial/test_cell_update_protocol, 0, 0, 0, 0, 16, 0
   spatial/test_diamond, 0, 0, 0, 0, 13, 0
   spatial/test_pole_angular_closure, 0, 0, 0, 0, 28, 0
   test_convergence, 0, 0, 1, 0, 0, 0
   test_pending_ports, 5, 0, 0, 0, 0, 0
   test_vv_harness_audit, 9, 0, 0, 0, 0, 0

Equation coverage
-----------------

Every Sphinx ``.. math:: :label:`` block declared in ``docs/theory/*.rst`` and the number of tests carrying ``@pytest.mark.verifies("label")`` that reference it.

.. csv-table::
   :header: Equation label, Tests
   :widths: 50, 10

   ``matrix-eigenvalue``, 178
   ``mg-balance``, 168
   ``peierls-unified``, 158
   ``one-group-kinf``, 132
   ``peierls-rank-n-bc-closure``, 124
   ``reflective-bc``, 110
   ``alpha-recursion``, 100
   ``wdd-closure``, 100
   ``wdd-face``, 100
   ``collision-rate``, 94
   ``alpha-cylindrical``, 74
   ``mm-weights``, 74
   ``multigroup``, 67
   ``ki3-def``, 64
   ``e3-def``, 61
   ``self-slab``, 55
   ``self-cyl``, 54
   ``p-inf``, 53
   ``self-sph``, 52
   ``balance-general``, 51
   ``chord-length``, 51
   ``flux-moments``, 49
   ``wigner-seitz``, 49
   ``attenuation``, 48
   ``optical-thickness``, 48
   ``scalar-flux-integral``, 48
   ``cp-kernel-differential-identities``, 36
   ``flat-source``, 35
   ``cp-keff-update``, 34
   ``first-flight-kernel``, 34
   ``matrix-A-def``, 34
   ``matrix-B-def``, 34
   ``neutron-balance``, 34
   ``optical-path``, 34
   ``pcell-from-smat``, 34
   ``pin-from-reciprocity``, 34
   ``rcp-from-double-antideriv``, 34
   ``s-integral``, 34
   ``self-double-integral``, 34
   ``surface-to-region``, 34
   ``surface-to-surface``, 34
   ``keff-mean``, 33
   ``peierls-equation``, 33
   ``sigma-keff``, 33
   ``free-flight``, 31
   ``dd-slab``, 30
   ``transport-spherical``, 29
   ``chi-sampling``, 28
   ``decompose``, 28
   ``scattering-cdf``, 28
   ``transport-cylindrical``, 28
   ``cp-flat-source-derivation``, 27
   ``cp-flat-source-double-integral``, 27
   ``cp-unified-outer-integration``, 27
   ``fission-matrix``, 26
   ``removal-matrix``, 26
   ``azimuthal-angles``, 24
   ``dc-slab``, 24
   ``peierls-specular-bc-defn``, 24
   ``second-diff-cyl``, 24
   ``second-diff-sph``, 24
   ``effective-spacing``, 22
   ``peierls-cyl-3d-mode-formula``, 22
   ``pitch-recovery``, 22
   ``ray-circle``, 22
   ``en-kernel-derivative``, 20
   ``kin-kernel-derivative``, 20
   ``peierls-rank-n-stability``, 20
   ``dd-cartesian-1d``, 17
   ``direction-sampling``, 16
   ``discrete-measure-integrate``, 16
   ``fission-weight``, 16
   ``keff-cycle``, 16
   ``roulette-conservation``, 16
   ``roulette-prob``, 16
   ``transport-cartesian``, 16
   ``peierls-greens-hollow-sph-architecture``, 15
   ``second-diff-general``, 15
   ``peierls-greens-annulus-architecture``, 14
   ``pn-scatter``, 14
   ``peierls-greens-slab-asym-architecture``, 13
   ``complementarity``, 12
   ``kinf-1g``, 12
   ``kinf-mg``, 12
   ``periodic-bc``, 12
   ``reciprocity``, 12
   ``ws-pitch``, 12
   ``inf-hom-balance``, 11
   ``two-group-A``, 11
   ``two-group-Ainv``, 11
   ``two-group-F``, 11
   ``two-group-M``, 11
   ``atalay-eq42-extrapolated-endpoint``, 10
   ``kll-1974-slab-flux``, 10
   ``kll-1974-sphere-flux``, 10
   ``peierls-greens-cylinder-architecture``, 10
   ``peierls-greens-slab-architecture``, 10
   ``peierls-vacuum-bc-flux``, 10
   ``peierls-vacuum-bc-row-sum-gate``, 10
   ``peierls-white-bc``, 10
   ``singular-eigenfunction-eq42``, 10
   ``dd-recurrence``, 9
   ``singular-eigenfunction-eq46``, 9
   ``tau-m``, 9
   ``tau-p``, 9
   ``bare-slab-buckling``, 8
   ``bare-slab-critical-equation``, 8
   ``bare-slab-eigenfunction``, 8
   ``cp-inner-integral-antiderivative``, 8
   ``diffusion-M-matrix``, 8
   ``diffusion-back-substitution``, 8
   ``diffusion-coefficient``, 8
   ``diffusion-exponential-branch``, 8
   ``diffusion-interface-matching``, 8
   ``diffusion-matching-matrix``, 8
   ``diffusion-mode-decomposition``, 8
   ``diffusion-operator``, 8
   ``diffusion-region-ode``, 8
   ``diffusion-spurious-root-validation``, 8
   ``diffusion-transcendental``, 8
   ``diffusion-trigonometric-branch``, 8
   ``dd-solve``, 7
   ``gauss-legendre-visibility-cone``, 7
   ``bar-psi``, 6
   ``boyd-eq-45``, 6
   ``characteristic-ode``, 6
   ``kin-kernel-special-values``, 6
   ``peierls-greens-slab-T``, 6
   ``sn-streaming-reciprocity``, 6
   ``cp-outer-integral-antiderivative``, 5
   ``dd-cartesian-2d``, 5
   ``en-kernel-special-values``, 5
   ``peierls-cyl-Gbc-3d-final``, 5
   ``peierls-greens-cylinder-T``, 5
   ``singular-eigenfunction-eq40``, 5
   ``transport-cartesian-2d``, 5
   ``xs-interp``, 5
   ``absorption-xs``, 4
   ``atalay-eq46-slab-eq54-sphere-parity``, 4
   ``en-kernel-integral``, 4
   ``fission-source``, 4
   ``fixed-source-solve``, 4
   ``hebert-3-323``, 4
   ``keff-update``, 4
   ``macro-sum``, 4
   ``peierls-greens-V-alpha-1``, 4
   ``peierls-greens-V-alpha-2``, 4
   ``peierls-vacuum-bc-slab``, 4
   ``peierls-white-bc-slab``, 4
   ``sn-mms-hetero-psi``, 4
   ``sn-mms-hetero-qext``, 4
   ``two-group-charpoly``, 4
   ``two-group-roots``, 4
   ``atalay-table2-slab-reflected-isotropic``, 3
   ``atalay-table2-slab-reflected-r099-precision-floor``, 3
   ``atalay-table2-slab-vacuum-isotropic``, 3
   ``atalay-table3-slab-vacuum-anisotropic``, 3
   ``cp-escape-from-p-cell``, 3
   ``dd-curvilinear-scalar``, 3
   ``dd-mm-closure-constants``, 3
   ``dd-slab-scalar``, 3
   ``delta-psi``, 3
   ``hetero-tolerance``, 3
   ``isotropic-source``, 3
   ``moc-keff-update``, 3
   ``moc-mms-psi-ref``, 3
   ``moc-mms-qext``, 3
   ``moc-wigner-seitz``, 3
   ``number-density``, 3
   ``peierls-vacuum-bc-cylinder``, 3
   ``peierls-vacuum-bc-sphere``, 3
   ``sigma-zero``, 3
   ``atalay-eq54-sphere-vacuum-isotropic``, 2
   ``atalay-table6-eigenvalue-moderate-d-consistency``, 2
   ``cp-second-difference-operator``, 2
   ``peierls-greens-slab-asym-method-of-images``, 2
   ``peierls-slab-Gbc-mode``, 2
   ``peierls-slab-Pesc-mode``, 2
   ``richardson-diffusion``, 2
   ``roulette-restore``, 2
   ``singular-eigenfunction-eq54``, 2
   ``sn-case-back-substitution``, 2
   ``sn-case-matching-matrix``, 2
   ``sn-case-per-ordinate``, 2
   ``sn-case-physical-validation``, 2
   ``sn-case-real-basis``, 2
   ``sn-case-slope-matrix``, 2
   ``sn-case-spatial-modes``, 2
   ``sn-mms-2d-2g-psi``, 2
   ``sn-mms-cylindrical-aniso-psi``, 2
   ``sn-mms-cylindrical-aniso-qext``, 2
   ``sn-mms-p1-qext``, 2
   ``sn-mms-spherical-aniso-psi``, 2
   ``sn-mms-spherical-aniso-qext``, 2
   ``addition-theorem``, 1
   ``branching``, 1
   ``collision-estimator``, 1
   ``e1-decomposition``, 1
   ``majorant``, 1
   ``nm-1980-reflected-slab-fn``, 1
   ``normalisation``, 1
   ``peierls-class-b-Pss-homogeneous``, 1
   ``peierls-greens-T00-integrand``, 1
   ``peierls-greens-V-alpha-3``, 1
   ``peierls-greens-annulus-3d-chord-scaling``, 1
   ``peierls-greens-annulus-impact-parameter-partition``, 1
   ``peierls-greens-annulus-through-rank2``, 1
   ``peierls-greens-cylinder-bounce-period``, 1
   ``peierls-greens-cylinder-impact-parameter``, 1
   ``peierls-greens-cylinder-in-plane-speed``, 1
   ``peierls-greens-cylinder-trajectory``, 1
   ``peierls-greens-function-architecture``, 1
   ``peierls-greens-hollow-sph-impact-parameter-partition``, 1
   ``peierls-greens-hollow-sph-through-rank2``, 1
   ``peierls-greens-slab-asym-closure``, 1
   ``peierls-greens-slab-asym-resolvent``, 1
   ``peierls-greens-slab-trajectory``, 1
   ``peierls-greens-surface-fixed-point``, 1
   ``real-spherical-harmonics``, 1
   ``sigT-computed``, 1
   ``singular-eigenfunction-eq5``, 1
   ``sn-mms-2d-2g-qext``, 1
   ``sn-mms-2d-psi``, 1
   ``sn-mms-2d-qext``, 1
   ``sn-mms-cylindrical-psi``, 1
   ``sn-mms-cylindrical-qext``, 1
   ``sn-mms-p1-psi``, 1
   ``sn-mms-psi``, 1
   ``sn-mms-qext``, 1
   ``sn-mms-spherical-psi``, 1
   ``sn-mms-spherical-qext``, 1
   ``splitting``, 1

Orphan equations
----------------

Equations with zero tests carrying ``@pytest.mark.verifies("label")``, excluding labels explicitly marked ``:vv-status: documented``. **29** of the testable equations found on theory pages are orphan.

- ``affine-bc-form``
- ``bc-tensor-decomposition``
- ``billiard-rank2-S``
- ``billiard-rank2-T``
- ``billiard-reflection-law``
- ``billiard-resolvent-neumann``
- ``billiard-transfer-operator``
- ``billiard-variant-alpha-rank1``
- ``fn-method-moment-space-AB-defs``
- ``fn-method-moment-space-bc-vacuum``
- ``fn-method-moment-space-fn-ansatz``
- ``fn-method-moment-space-galerkin-system``
- ``galerkin-spectral-bte``
- ``galerkin-spectral-carlvik-integral``
- ``galerkin-spectral-eq4``
- ``galerkin-spectral-orthogonality``
- ``inflow-mask-discrete``
- ``ordinate-partition-inflow-outflow``
- ``pole-mm-recurrence``
- ``spectrum-case-eigenfunction-equation``
- ``spectrum-case-eigenfunction-explicit``
- ``spectrum-continuum-eigenfunction``
- ``spectrum-dispersion-relation``
- ``spectrum-expansion-theorem``
- ``spectrum-full-decomposition``
- ``spectrum-transport-equation``
- ``spectrum-x-function``
- ``trace-sign-predicate``
- ``vacuum-legacy-vs-trace-correct``

Documented-only equations
-------------------------

Theory labels marked ``.. vv-status: <label> documented`` in their RST source. These are excluded from the orphan-equation gate because they are either definitional (no single implementing function — e.g. ``boltzmann``), describe a module whose Python port does not yet exist (e.g. the thermal-hydraulics / fuel-behaviour / reactor-kinetics equations), or have a deliberately deferred test paired with a tracking issue. **183** labels carry the directive. See ``docs/testing/architecture.rst``:ref:`vv-status-documented` for the full taxonomy.

- ``bailey-dome-recursion``
- ``bc-rank-n-as-sum-of-products``
- ``bessel-wronskian``
- ``bickley-integral``
- ``boltzmann``
- ``bundle-measure-disintegration``
- ``burst-criterion``
- ``c-in-remapping``
- ``case-dispersion-function``
- ``clad-heat``
- ``conservative-form``
- ``convergence-rate``
- ``coolant-energy``
- ``coolant-feedback``
- ``coolant-rate``
- ``creep-rate``
- ``dd-2d-balance-form``
- ``diagonal-operator-action``
- ``discrete-measure-definition``
- ``discrete-measure-partition``
- ``discrete-measure-pushforward``
- ``doppler-feedback``
- ``e1-small-tau-expansion``
- ``fb-bc4-displacement``
- ``fb-clad-strain``
- ``fb-fuel-heat``
- ``fb-fuel-strain``
- ``fb-swelling``
- ``fn-Fk-closed-forms``
- ``fn-Fk-integration-by-parts``
- ``fn-product-simpson-weights``
- ``fn-slab-B-long-division``
- ``fn-slab-collocation``
- ``fn-unified-matrix-entry``
- ``fn-x-function``
- ``fuel-heat``
- ``fuel-rate``
- ``galerkin-construction``
- ``galerkin-pair``
- ``galerkin-self-adjoint``
- ``gap-closure-event``
- ``gap-conductance``
- ``gas-pressure``
- ``group-flux``
- ``group-xs``
- ``kll-1974-slab-phi``
- ``kll-1974-sphere-phi``
- ``maxwellian``
- ``mode-conservation-target``
- ``morel-montry-clamp``
- ``nm1980-eq16-tau-zero``
- ``octant-direct-sum-tensor-product``
- ``octant-sign-predicate``
- ``one-over-E``
- ``operator-apply``
- ``operator-apply-transpose``
- ``operator-eigenvalue``
- ``operator-fixed-source``
- ``operator-solve``
- ``partition-round-trip``
- ``peierls-3d``
- ``peierls-M-rank-1``
- ``peierls-M-rank-2``
- ``peierls-WM-WL-asymmetric``
- ``peierls-bc-general``
- ``peierls-bc-operator``
- ``peierls-boltzmann``
- ``peierls-change-of-basis``
- ``peierls-class-b-Jn-canonical``
- ``peierls-class-b-hebert-closure``
- ``peierls-cyl-Pss-derivation``
- ``peierls-cyl-foundations``
- ``peierls-cylinder-equation``
- ``peierls-cylinder-green-2d``
- ``peierls-cylinder-nystrom``
- ``peierls-cylinder-polar``
- ``peierls-cylinder-r-prime``
- ``peierls-cylinder-ray-optical-depth``
- ``peierls-cylinder-rho-max``
- ``peierls-cylinder-row-sum-identity``
- ``peierls-davison-urho``
- ``peierls-delta-tracking-equivalence``
- ``peierls-e1-derivation``
- ``peierls-exp-stretched-mu``
- ``peierls-factored-kernel``
- ``peierls-greens-A1-split``
- ``peierls-greens-A5-specular``
- ``peierls-greens-L0``
- ``peierls-greens-Lp``
- ``peierls-greens-T-alpha``
- ``peierls-greens-T-mu-surf``
- ``peierls-greens-bounce-period-integral``
- ``peierls-greens-bounce-sum-alpha``
- ``peierls-greens-defining-bvp``
- ``peierls-greens-fixed-source-iteration``
- ``peierls-greens-garcia-convention``
- ``peierls-greens-hollow-sph-outer-only-resolvent``
- ``peierls-greens-k-inf``
- ``peierls-greens-mg-source``
- ``peierls-greens-mr-piecewise-tau``
- ``peierls-greens-mr-trajectory-segments``
- ``peierls-greens-mu-surf``
- ``peierls-greens-sanchez-A6``
- ``peierls-greens-slab-asym-monodromy``
- ``peierls-greens-slab-bounce-period``
- ``peierls-greens-trajectory-integral``
- ``peierls-greens-unification-resolvent``
- ``peierls-half-range-inner-products``
- ``peierls-integral-form``
- ``peierls-kernel-decomposition``
- ``peierls-ki1-derivation``
- ``peierls-mg-operator``
- ``peierls-operator-factorisation``
- ``peierls-operator-form``
- ``peierls-point-kernel-3d``
- ``peierls-polar-jacobian-cancellation``
- ``peierls-rank-n-P-esc-moment``
- ``peierls-rank-n-jacobian-derivation``
- ``peierls-scaled-chebyshev``
- ``peierls-slab-bare-critical``
- ``peierls-slab-foundations``
- ``peierls-slab-polar``
- ``peierls-specular-M-tridiagonal``
- ``peierls-specular-R-formula``
- ``peierls-specular-T-matrix``
- ``peierls-specular-multibounce-formula``
- ``peierls-sph-ps1982-foundations``
- ``peierls-sphere-G-bc``
- ``peierls-sphere-equation``
- ``peierls-sphere-green-3d``
- ``peierls-sphere-nystrom``
- ``peierls-sphere-polar``
- ``peierls-sphere-r-prime``
- ``peierls-sphere-ray-optical-depth``
- ``peierls-sphere-rho-max``
- ``peierls-sphere-row-sum-identity``
- ``peierls-surface-centred-chord-discriminant``
- ``peierls-surface-centred-chord-quadratic``
- ``peierls-surface-centred-tangent-angles``
- ``peierls-svd``
- ``peierls-tau-coordinate-transform``
- ``peierls-tensor-G-definition``
- ``peierls-tensor-P-definition``
- ``peierls-unified``
- ``per-face-inflow-mask``
- ``petrov-galerkin-construction``
- ``pi-r-equals-4pi-i``
- ``power-equation``
- ``precursor-equation``
- ``quadrature-selection-criterion``
- ``real-sh-addition-theorem``
- ``real-sh-discrete-orthogonality``
- ``real-sh-l0``
- ``real-sh-l1``
- ``real-sh-l2plus``
- ``scattering-as-tensor-product-sum``
- ``sigs-convention``
- ``singular-eigenfunction-eq5``
- ``sood-eq18-1g-balance``
- ``sood-eq19-kinf-1g``
- ``sood-eq20-kinf-1g-c-form``
- ``sood-eq25-2g-matrix``
- ``sood-eq29-kinf-2g-no-upscatter``
- ``sood-eq32-phi-ratio``
- ``sood-eq76-kinf-mg``
- ``streaming-as-tensor-product-sum``
- ``streaming-inverse-direct-sum``
- ``subgroup-of-o3-containment``
- ``sum-of-tensor-products``
- ``tensor-product-action``
- ``tensor-product-adjoint-distributivity``
- ``tensor-product-axis-wise-composition``
- ``tensor-product-inverse``
- ``trace-half-decomposition``
- ``transport-equation``
- ``vacuum-bc``
- ``wm72-coupled-linear-system``
- ``wm72-eq30-bare``
- ``wm72-eq31``
- ``wm72-eq32``
- ``wm72-q-formula``
- ``wm72-rho-bare-cylinder``
- ``wm72-singular-subtraction``

L0 error-catalog coverage
-------------------------

Every ``ERR-NNN`` entry in ``.claude/skills/vv-principles/error_catalog.md`` and the tests that carry ``@pytest.mark.catches("ERR-NNN")`` to guard it. A missing catcher is a publication-blocker for the error catalog.

.. csv-table::
   :header: Error tag, Catching tests
   :widths: 15, 10

   ``ERR-001``, 1
   ``ERR-002``, 1
   ``ERR-003``, 8
   ``ERR-004``, 1
   ``ERR-005``, 1
   ``ERR-006``, 2
   ``ERR-007``, 3
   ``ERR-008``, 1
   ``ERR-009``, 9
   ``ERR-010``, 1
   ``ERR-011``, 1
   ``ERR-012``, 1
   ``ERR-013``, 1
   ``ERR-014``, 1
   ``ERR-015``, 1
   ``ERR-016``, 2
   ``ERR-017``, 3
   ``ERR-018``, 1
   ``ERR-019``, 1
   ``ERR-020``, **0 (MISSING)**
   ``ERR-021``, 2
   ``ERR-022``, 1
   ``ERR-023``, 1
   ``ERR-024``, 1
   ``ERR-025``, 4
   ``ERR-026``, 4
   ``ERR-027``, 5
   ``ERR-028``, 1
   ``ERR-029``, 6
   ``ERR-030``, 2
   ``ERR-031``, **0 (MISSING)**
   ``ERR-032``, 4
   ``ERR-033``, 4
   ``ERR-034``, 1
   ``ERR-035``, 1
   ``ERR-036``, 8
   ``ERR-037``, 11
   ``ERR-038``, 5
   ``ERR-039``, 2
   ``ERR-040``, **0 (MISSING)**
   ``ERR-041``, **0 (MISSING)**
   ``ERR-042``, **0 (MISSING)**
   ``ERR-043``, 8
   ``ERR-044``, 4
   ``ERR-045``, **0 (MISSING)**
   ``ERR-046``, 8
   ``ERR-047``, **0 (MISSING)**

Unmarked tests
--------------

**28 tests** have no V&V level marker.
This is a gap — every test in the tree should carry either
a physics-ladder marker (``l0``..``l3``) or the orthogonal
``foundation`` marker (``@pytest.mark.foundation``) for
tests that verify software invariants rather than physics
equations. See ``docs/testing/architecture.rst``
:ref:`vv-foundation-tests` for the taxonomy.

.. csv-table::
   :header: File, Unmarked tests
   :widths: 60, 10

   ``tests/sn/regression/test_dd_regression.py``, 11
   ``tests/sn/test_boundary_conditions.py``, 11
   ``tests/derivations/test_atkinson_product_nystrom.py``, 3
   ``tests/derivations/test_fn_projection_vs_kll_flux.py``, 2
   ``tests/derivations/test_fn_la13511_slab_reflected.py``, 1

