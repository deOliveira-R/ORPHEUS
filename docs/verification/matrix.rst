Verification Matrix
===================

.. note::

   Auto-generated from ``tests._harness.registry.TEST_REGISTRY``
   by ``tools/verification/generate_matrix.py``. Do not edit by
   hand — changes will be overwritten on the next rebuild.

Total tests collected: **6202**

V&V level distribution
----------------------

.. csv-table::
   :header: Level, Count, Share
   :widths: 15, 10, 10

   L0, 1170, 18.9%
   L1, 1065, 17.2%
   L2, 49, 0.8%
   L3, 0, 0.0%
   foundation, 3908, 63.0%
   unmarked, 10, 0.2%

Tagging source
--------------

How each test acquired its V&V level (see ``tests/conftest.py`` for the precedence chain).

.. csv-table::
   :header: Source, Count
   :widths: 20, 10

   explicit, 6113
   verify, 0
   class-name, 46
   func-name, 0
   case, 33
   unmarked, 10

Module × level grid
-------------------

.. csv-table::
   :header: Module, L0, L1, L2, L3, FD, ??
   :widths: 40, 6, 6, 6, 6, 6, 6

   analytical/test_cp_standoff_curvilinear, 0, 2, 1, 0, 0, 0
   analytical/test_kinf_homogeneous, 0, 31, 0, 0, 0, 0
   analytical/test_kinf_homogeneous_tolerance, 0, 8, 0, 0, 0, 0
   analytical/test_l1_standoff_slab_cylinder, 0, 10, 0, 0, 0, 0
   analytical/test_mms_prescribed_inflow, 0, 3, 0, 0, 1, 0
   analytical/test_phase_c_crosscheck, 0, 8, 0, 0, 0, 0
   analytical/test_prescribed_inflow_consistency, 0, 0, 0, 0, 2, 0
   analytical/test_si_convergence_rate, 0, 7, 0, 0, 2, 0
   cartesian_2d/test_2d_full_field_oracle, 0, 0, 0, 0, 8, 0
   cartesian_2d/test_2d_l2_face_view_unit_source, 0, 0, 0, 0, 7, 0
   cartesian_2d/test_2d_l2_matvec_correctness, 0, 2, 0, 0, 2, 0
   cartesian_2d/test_2d_octant_sweep_equivalence, 0, 7, 0, 0, 0, 0
   cartesian_2d/test_discrete_ordinates_2d, 0, 0, 2, 0, 0, 0
   cartesian_2d/test_l2_boundary_face_view, 0, 0, 0, 0, 11, 0
   cartesian_2d/test_scan_march_equivalence, 0, 0, 0, 0, 11, 0
   core/test_affine_carve_baseline, 0, 0, 0, 0, 6, 0
   core/test_cell_balance_for_streaming, 0, 0, 0, 0, 9, 0
   core/test_cell_kernel_batch, 11, 0, 0, 0, 1, 2
   core/test_cell_visit_c_stamp, 0, 0, 0, 0, 3, 0
   core/test_dag_walk, 0, 0, 0, 0, 7, 0
   core/test_diamond, 0, 0, 0, 0, 53, 0
   core/test_discretization_scheme_protocol, 0, 0, 0, 0, 17, 0
   core/test_ordinate_scan, 52, 0, 0, 0, 0, 0
   core/test_ordinate_scan_joint_batch, 5, 0, 0, 0, 0, 0
   core/test_phase_c_gates, 8, 0, 0, 0, 9, 0
   core/test_reframe_moment_intent, 7, 0, 0, 0, 0, 0
   core/test_sweep_cache, 28, 0, 0, 0, 0, 0
   core/test_sweep_graph, 76, 0, 0, 0, 0, 0
   core/test_sweep_graph_nd_admission, 0, 0, 0, 0, 42, 0
   core/test_sweep_graph_window_equivalence, 0, 0, 0, 0, 20, 0
   core/test_sweep_regression, 10, 0, 0, 0, 0, 0
   core/test_sweep_schedule, 0, 0, 0, 0, 9, 0
   core/test_sweep_schedule_nd, 0, 0, 0, 0, 9, 0
   core/test_sweep_vs_apply_consistency, 0, 0, 0, 0, 57, 0
   core/test_transport_sweep_ng2_layout_guard, 0, 0, 0, 0, 2, 0
   core/test_unified_sweep_dispatch, 0, 0, 0, 0, 28, 0
   core/test_wavefront_cumprod_equivalence, 0, 0, 0, 0, 4, 0
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
   curvilinear/test_apply_matvec_cylinder_invariants, 24, 0, 0, 0, 0, 0
   curvilinear/test_compute_psi_half_per_level, 26, 0, 0, 0, 2, 0
   curvilinear/test_coupled_pole_mu_level_invariant, 0, 0, 0, 0, 12, 0
   curvilinear/test_cyl_sweep_regression, 4, 0, 4, 0, 0, 0
   curvilinear/test_pole_angular_closure, 0, 0, 0, 0, 13, 0
   curvilinear/test_psi_half_angle_seed, 4, 1, 0, 0, 20, 0
   curvilinear/test_si_cyl_20cell_nan_regression, 0, 0, 0, 0, 4, 0
   curvilinear/test_sph_sweep_regression, 13, 0, 0, 0, 0, 0
   curvilinear/test_streaming_equilibrium_curvilinear, 27, 0, 0, 0, 0, 0
   curvilinear/test_tau_producer_equivalence, 0, 0, 0, 0, 5, 0
   curvilinear/test_unified_matvec_cylinder, 29, 2, 0, 0, 0, 0
   curvilinear/test_unified_matvec_sphere, 2, 0, 0, 0, 0, 0
   curvilinear/test_w1_clamp_silent_on_flat, 0, 2, 0, 0, 2, 0
   data/test_chi_invariant_enforcement, 0, 0, 0, 0, 13, 0
   data/test_chi_mix_production_weighting, 0, 0, 0, 0, 8, 0
   data/test_cross_section_data, 11, 0, 0, 0, 0, 0
   data/test_emission_spectrum, 0, 0, 0, 0, 15, 0
   data/test_energy_grid, 0, 0, 0, 0, 28, 0
   data/test_gendf_canonical_order, 0, 0, 0, 0, 6, 0
   data/test_group_permutation_invariance, 0, 0, 10, 0, 0, 0
   data/test_mixture, 4, 0, 0, 0, 0, 0
   data/test_mixture_condense, 0, 0, 0, 0, 38, 0
   data/test_mixture_scattering_ratio, 0, 0, 0, 0, 3, 0
   data/test_mixture_transport_xs, 0, 1, 0, 0, 4, 0
   data/test_mixture_xs_balance, 0, 0, 0, 0, 75, 0
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
   derivations/test_continuous_registry_lazy, 0, 0, 0, 0, 5, 0
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
   derivations/test_peierls_fission_source_indexing, 0, 3, 0, 0, 0, 0
   derivations/test_peierls_geometry, 0, 0, 0, 0, 32, 0
   derivations/test_peierls_greens_function_annulus_solver, 0, 13, 0, 0, 0, 0
   derivations/test_peierls_greens_function_annulus_symbolic, 0, 0, 0, 0, 22, 0
   derivations/test_peierls_greens_function_cylinder_mr, 0, 5, 0, 0, 5, 0
   derivations/test_peierls_greens_function_cylinder_mr_xverif, 0, 1, 0, 0, 0, 0
   derivations/test_peierls_greens_function_cylinder_solver, 0, 13, 0, 0, 0, 0
   derivations/test_peierls_greens_function_cylinder_symbolic, 0, 0, 0, 0, 12, 0
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
   derivations/test_sn_mms_ld_2d_stress_symbolic, 0, 0, 0, 0, 8, 0
   derivations/test_sn_mms_nonvacuum_symbolic, 0, 0, 0, 0, 9, 0
   derivations/test_sood_registry_cache, 0, 0, 0, 0, 15, 0
   derivations/test_sood_registry_compatibility, 0, 2, 0, 0, 107, 0
   derivations/test_sood_registry_wide_bare_critical, 0, 17, 0, 0, 2, 0
   derivations/test_sood_registry_wide_kinf, 0, 0, 0, 0, 49, 0
   derivations/test_trajectory_resolvent_billiard, 0, 0, 0, 0, 11, 0
   derivations/test_trajectory_resolvent_chord_oracle, 0, 0, 0, 0, 21, 0
   derivations/test_trajectory_resolvent_power_iterate, 0, 0, 0, 0, 6, 0
   derivations/test_xs_library_validation, 0, 0, 0, 0, 2, 0
   diffusion/test_augmented_mesh, 0, 0, 0, 0, 12, 0
   diffusion/test_boundary_realizer, 4, 0, 0, 0, 27, 0
   diffusion/test_continuous_reference, 0, 7, 0, 0, 0, 0
   diffusion/test_mms, 0, 2, 0, 0, 0, 0
   diffusion/test_operators, 8, 0, 0, 0, 22, 0
   diffusion/test_properties, 3, 0, 0, 0, 0, 0
   diffusion/test_solver, 0, 0, 3, 0, 15, 0
   eigenvalue/test_heterogeneous_transport, 0, 2, 0, 0, 0, 0
   eigenvalue/test_keff_2d, 19, 0, 0, 0, 0, 0
   eigenvalue/test_keff_curvilinear, 0, 19, 12, 0, 0, 0
   eigenvalue/test_keff_estimator_gate, 0, 0, 0, 0, 7, 0
   eigenvalue/test_keff_slab, 0, 6, 0, 0, 0, 0
   fields/test_angular_boundary_flux, 0, 0, 0, 0, 36, 0
   fields/test_angular_boundary_source_sink_residual, 0, 0, 0, 0, 30, 0
   fields/test_angular_flux, 0, 0, 0, 0, 24, 0
   fields/test_coefficient_fields, 0, 0, 0, 0, 9, 0
   fields/test_scalar_boundary_flux, 0, 0, 0, 0, 15, 0
   frames/test_harmonic_frame, 0, 0, 0, 0, 11, 0
   geometry/test_bc_equivalence_snapshot, 0, 8, 0, 0, 0, 0
   geometry/test_bc_errors, 0, 0, 0, 0, 11, 0
   geometry/test_bc_universal_invariants, 0, 30, 0, 0, 11, 0
   geometry/test_bound_compat, 10, 0, 0, 0, 0, 0
   geometry/test_boundary, 0, 0, 0, 0, 25, 0
   geometry/test_boundary_trace_law, 0, 0, 0, 0, 15, 0
   geometry/test_geometry, 0, 0, 0, 0, 55, 0
   geometry/test_law_composition, 0, 2, 0, 0, 16, 0
   geometry/test_mesh, 0, 0, 0, 0, 10, 0
   geometry/test_reduced_operator, 0, 0, 0, 0, 47, 0
   geometry/test_structured_geometry, 0, 0, 0, 0, 36, 0
   homogeneous/test_continuous_reference, 0, 9, 0, 0, 0, 0
   homogeneous/test_homogeneous, 0, 14, 0, 0, 0, 0
   mc/test_convergence, 0, 0, 3, 0, 0, 0
   mc/test_cross_verification, 0, 0, 2, 0, 0, 0
   mc/test_gaps, 7, 9, 0, 0, 0, 0
   mc/test_monte_carlo, 0, 12, 0, 0, 0, 0
   mc/test_properties, 24, 0, 0, 0, 0, 0
   mms/test_curvilinear_aniso_convergence, 0, 7, 0, 0, 0, 0
   mms/test_curvilinear_aniso_scattering_p1, 2, 0, 0, 0, 0, 0
   mms/test_curvilinear_operator_admits_anisotropic_mms, 0, 2, 0, 0, 0, 0
   mms/test_curvilinear_operator_admits_mms, 0, 2, 0, 0, 0, 0
   mms/test_curvilinear_pole_cell_characterization, 0, 4, 0, 0, 0, 0
   mms/test_ld_2d_boundary_promise, 0, 10, 0, 0, 1, 0
   mms/test_mms, 0, 2, 0, 0, 0, 0
   mms/test_mms_2d, 0, 3, 0, 0, 0, 0
   mms/test_mms_aniso, 0, 2, 0, 0, 0, 0
   mms/test_mms_curvilinear, 0, 2, 0, 0, 0, 0
   mms/test_mms_heterogeneous, 0, 4, 0, 0, 0, 0
   mms/test_mms_ld_2d, 0, 11, 0, 0, 12, 1
   mms/test_mms_ld_slab, 0, 5, 0, 0, 2, 0
   mms/test_space_angle_separability, 0, 6, 0, 0, 0, 0
   moc/test_mms, 0, 3, 0, 0, 0, 0
   moc/test_moc, 0, 3, 0, 0, 0, 0
   moc/test_properties, 4, 0, 0, 0, 0, 0
   moc/test_quadrature, 24, 0, 0, 0, 0, 0
   moc/test_ray_tracing, 22, 0, 0, 0, 0, 0
   moc/test_verification, 27, 15, 6, 0, 0, 0
   numerics/test_affine_flux_algebra, 0, 0, 0, 0, 34, 0
   numerics/test_angular_trace_space, 10, 5, 0, 0, 12, 0
   numerics/test_diagonal_operator, 19, 0, 0, 0, 3, 0
   numerics/test_eigenvalue, 0, 39, 0, 0, 0, 0
   numerics/test_estimators_as_functionals, 0, 0, 0, 0, 4, 0
   numerics/test_face_layout, 0, 0, 0, 0, 15, 0
   numerics/test_field, 0, 0, 0, 0, 22, 0
   numerics/test_frame, 0, 0, 0, 0, 18, 0
   numerics/test_full_field_space, 0, 0, 0, 0, 6, 0
   numerics/test_green_operator, 0, 0, 0, 0, 11, 0
   numerics/test_incoming_ordinate_mask_tensor, 13, 0, 0, 0, 0, 0
   numerics/test_indicator_basis, 0, 0, 0, 0, 11, 0
   numerics/test_inner_product_functional, 4, 0, 0, 0, 0, 0
   numerics/test_inverse_universal, 0, 0, 0, 0, 24, 0
   numerics/test_iteration, 0, 1, 0, 0, 18, 0
   numerics/test_matrix_inverse_operator, 0, 0, 0, 0, 24, 0
   numerics/test_measure, 0, 16, 0, 0, 32, 0
   numerics/test_measure_partition, 12, 0, 0, 0, 0, 0
   numerics/test_measure_phase, 0, 0, 0, 0, 11, 0
   numerics/test_operator, 0, 0, 0, 0, 58, 0
   numerics/test_operator_capability_predicates, 0, 0, 0, 0, 27, 0
   numerics/test_operator_protocols, 0, 0, 0, 0, 14, 0
   numerics/test_outer_dyad, 9, 0, 0, 0, 0, 0
   numerics/test_periodic_wrap_operator, 5, 0, 0, 0, 0, 0
   numerics/test_permutation_operator, 9, 2, 0, 0, 0, 0
   numerics/test_quadrature_directional, 0, 0, 0, 0, 23, 0
   numerics/test_registry, 0, 0, 0, 0, 37, 0
   numerics/test_registry_mixin, 0, 0, 0, 0, 10, 0
   numerics/test_rules_1d, 0, 5, 0, 0, 16, 0
   numerics/test_rules_product, 0, 3, 0, 0, 14, 0
   numerics/test_rules_sphere, 0, 7, 0, 0, 29, 0
   numerics/test_space, 0, 0, 0, 0, 15, 0
   numerics/test_space_algebra, 0, 0, 0, 0, 18, 0
   numerics/test_spatial_moment_space, 0, 0, 0, 0, 22, 0
   numerics/test_spherical_harmonic_basis, 4, 7, 0, 0, 0, 0
   numerics/test_spherical_harmonic_space, 0, 8, 0, 0, 6, 0
   numerics/test_symmetry, 0, 0, 0, 0, 71, 0
   numerics/test_tensor_product_operator, 25, 0, 0, 0, 0, 0
   numerics/test_vector_protocol, 0, 0, 0, 0, 8, 0
   numerics/test_weighted_indicator_basis, 0, 0, 0, 0, 9, 0
   operators/test_angular_average_operator, 12, 4, 0, 0, 0, 0
   operators/test_apply_full_field_codomain, 0, 0, 0, 0, 15, 0
   operators/test_bc_extraction_2d, 2, 3, 0, 0, 3, 0
   operators/test_bc_extraction_matvec, 3, 0, 0, 0, 30, 0
   operators/test_boundary_conditions, 0, 0, 0, 0, 11, 0
   operators/test_capability_survival, 0, 0, 0, 0, 10, 0
   operators/test_collision_operator, 0, 0, 0, 0, 54, 0
   operators/test_fission_adjoint, 0, 0, 0, 0, 8, 0
   operators/test_fission_kernel_crosscheck, 0, 0, 0, 0, 6, 0
   operators/test_fission_operator, 0, 0, 0, 0, 18, 0
   operators/test_frame_conjugate_carve, 0, 0, 0, 0, 11, 0
   operators/test_g_adjoint_reciprocity, 0, 0, 0, 0, 12, 0
   operators/test_green_operator_sn, 0, 0, 0, 0, 5, 0
   operators/test_inverse_operator_equivalence, 0, 0, 0, 0, 7, 0
   operators/test_invertible_operator, 1, 10, 0, 0, 22, 0
   operators/test_isotropic_scattering, 0, 0, 0, 0, 16, 0
   operators/test_legendre_moment_scattering, 9, 0, 0, 0, 0, 0
   operators/test_loss_action_convention, 0, 0, 0, 0, 4, 0
   operators/test_native_matvec, 0, 0, 0, 0, 18, 0
   operators/test_one_octant_walk, 0, 0, 0, 0, 3, 0
   operators/test_one_representation_instance, 0, 0, 0, 0, 2, 0
   operators/test_operator_block_role, 0, 0, 0, 0, 20, 0
   operators/test_operators_apply_typed, 0, 0, 0, 0, 17, 0
   operators/test_pure_L_sigma_free, 0, 0, 0, 0, 9, 0
   operators/test_removal_form_matvec_sweep, 0, 0, 0, 0, 20, 0
   operators/test_scattering_adjoint, 0, 0, 0, 0, 15, 0
   operators/test_scattering_kernel_crosscheck, 0, 0, 0, 0, 4, 0
   operators/test_scattering_operator, 1, 0, 0, 0, 69, 0
   operators/test_sn_boundary_operator, 0, 0, 0, 0, 20, 0
   operators/test_sn_boundary_realizer, 0, 20, 0, 0, 0, 0
   operators/test_snmesh_realizer_wiring, 0, 11, 0, 0, 0, 0
   operators/test_solver_components, 27, 0, 0, 0, 0, 0
   operators/test_streaming_operator, 0, 0, 0, 0, 54, 0
   operators/test_streaming_operator_decomposition, 21, 0, 0, 0, 0, 0
   operators/test_typed_residual_evaluation, 1, 0, 0, 0, 6, 0
   primitives/test_axis_native_construction, 0, 0, 0, 0, 15, 0
   primitives/test_axis_primitive, 0, 0, 0, 0, 23, 0
   primitives/test_boundary_face_layout, 0, 0, 0, 0, 5, 0
   primitives/test_cell_flattening_invariant, 0, 0, 0, 0, 3, 0
   primitives/test_dag_ownership, 0, 0, 0, 0, 18, 0
   primitives/test_face_name_crosswalk, 0, 0, 0, 0, 4, 0
   primitives/test_harmonic_moment_flux, 0, 0, 0, 0, 31, 0
   primitives/test_method_space, 5, 0, 0, 0, 0, 0
   primitives/test_octants_property, 60, 0, 0, 0, 0, 0
   primitives/test_properties, 4, 0, 0, 0, 0, 0
   primitives/test_quadrature, 49, 0, 0, 0, 0, 0
   primitives/test_snmesh_consumes_reduced, 0, 0, 0, 0, 16, 0
   primitives/test_snmesh_materials_pr_typed_0, 0, 0, 0, 0, 7, 0
   primitives/test_solution, 0, 0, 0, 0, 31, 0
   primitives/test_typed_source_sinks, 0, 0, 0, 0, 37, 0
   regression/test_dd_regression, 0, 0, 0, 0, 13, 0
   residuals/test_typed_residuals, 0, 0, 0, 0, 33, 0
   slab/test_dd_recurrence, 1, 0, 0, 0, 0, 0
   slab/test_unified_matvec_slab, 2, 2, 0, 0, 0, 0
   sn/test_condensation, 0, 10, 0, 0, 0, 0
   sn/test_homogenization, 12, 0, 0, 0, 0, 0
   sn/test_material_xs_field_typed, 0, 0, 0, 0, 10, 0
   solve/test_2d_anisotropic_windowing, 0, 6, 0, 0, 0, 0
   solve/test_affine_carve_bit_identity, 0, 0, 0, 0, 3, 0
   solve/test_b1pp_verification, 6, 3, 0, 0, 0, 0
   solve/test_d3_admission, 0, 5, 0, 0, 2, 0
   solve/test_fixed_source_2d_equivalence, 0, 2, 0, 0, 0, 0
   solve/test_fixed_source_g1, 0, 5, 0, 0, 0, 0
   solve/test_flux_displacement_diagnostics, 0, 4, 0, 0, 0, 0
   solve/test_gauss_seidel_reification, 0, 0, 0, 0, 7, 0
   solve/test_krylov_curvilinear_precond_safety, 0, 4, 0, 0, 0, 0
   solve/test_krylov_restart_signature, 0, 12, 0, 0, 0, 0
   solve/test_scan_march_end_to_end, 0, 4, 0, 0, 0, 0
   solve/test_seed_threading_spy, 0, 0, 0, 0, 1, 0
   solve/test_si_gate_dispatch, 0, 0, 0, 0, 4, 0
   solve/test_si_single_primitive_contract, 0, 0, 0, 0, 2, 0
   spatial/test_affine_closure, 0, 0, 0, 0, 3, 0
   spatial/test_ld_slope_frame, 0, 1, 0, 0, 1, 0
   spatial/test_ld_ubld_primitive, 0, 0, 0, 0, 11, 0
   spatial/test_ld_ubld_symbolic, 0, 0, 0, 0, 6, 0
   spatial/test_linear_discontinuous, 0, 0, 0, 0, 20, 1
   spatial/test_moment_axis_predicates, 0, 0, 0, 0, 6, 0
   spatial/test_ordinate_scan_reset, 3, 2, 0, 0, 0, 0
   spatial/test_pairing_diffusion_limit, 0, 0, 0, 0, 6, 0
   spatial/test_scheme_reaction_rate_contract, 0, 0, 0, 0, 10, 0
   spatial/test_spatial_moment_field_space, 0, 0, 0, 0, 12, 0
   test_convergence, 0, 0, 1, 0, 0, 0
   test_layer_imports, 0, 0, 0, 0, 310, 0
   test_pending_ports, 5, 0, 0, 0, 0, 0
   test_pyright_ratchet, 0, 0, 0, 0, 1, 0
   test_vv_harness_audit, 9, 0, 0, 0, 0, 0
   transport/test_field_units, 0, 0, 0, 0, 43, 0
   transport/test_full_field, 0, 0, 0, 0, 13, 0
   transport/test_functional_category, 0, 0, 0, 0, 11, 0
   transport/test_integral_kernel_category, 0, 0, 0, 0, 13, 0
   transport/test_integrated_reaction_rate, 0, 0, 0, 0, 6, 0
   transport/test_material_mesh, 0, 0, 0, 0, 12, 0
   transport/test_method, 0, 0, 0, 0, 4, 0
   transport/test_multiplication_operator, 0, 0, 0, 0, 22, 0
   transport/test_reaction_rate_functional, 0, 0, 0, 0, 7, 0
   transport/test_timed_full_field, 0, 0, 0, 0, 38, 0

Equation coverage
-----------------

Every Sphinx ``.. math:: :label:`` block declared in ``docs/theory/*.rst`` and the number of tests carrying ``@pytest.mark.verifies("label")`` that reference it.

.. csv-table::
   :header: Equation label, Tests
   :widths: 50, 10

   ``matrix-eigenvalue``, 227
   ``mg-balance``, 182
   ``peierls-unified``, 158
   ``one-group-kinf``, 146
   ``peierls-rank-n-bc-closure``, 124
   ``reflective-bc``, 112
   ``alpha-recursion``, 102
   ``wdd-closure``, 102
   ``wdd-face``, 102
   ``collision-rate``, 94
   ``multigroup``, 81
   ``alpha-cylindrical``, 75
   ``mm-weights``, 75
   ``ki3-def``, 64
   ``fission-matrix``, 63
   ``removal-matrix``, 63
   ``e3-def``, 61
   ``blelloch-1990-eq-1-5``, 57
   ``self-slab``, 55
   ``hebert-3-432``, 54
   ``self-cyl``, 54
   ``balance-general``, 53
   ``p-inf``, 53
   ``self-sph``, 52
   ``chord-length``, 51
   ``flux-moments``, 51
   ``transport-cartesian``, 51
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
   ``transport-spherical``, 33
   ``dd-slab``, 32
   ``free-flight``, 31
   ``transport-cylindrical``, 30
   ``chi-sampling``, 28
   ``decompose``, 28
   ``scattering-cdf``, 28
   ``sn-curvilinear-homogeneous-kinf-recovery``, 28
   ``cp-flat-source-derivation``, 27
   ``cp-flat-source-double-integral``, 27
   ``cp-unified-outer-integration``, 27
   ``azimuthal-angles``, 24
   ``dc-slab``, 24
   ``loss-rep-resolution-a``, 24
   ``peierls-specular-bc-defn``, 24
   ``second-diff-cyl``, 24
   ``second-diff-sph``, 24
   ``inf-hom-balance``, 23
   ``two-group-A``, 23
   ``two-group-Ainv``, 23
   ``two-group-F``, 23
   ``two-group-M``, 23
   ``effective-spacing``, 22
   ``peierls-cyl-3d-mode-formula``, 22
   ``pitch-recovery``, 22
   ``pn-scatter``, 22
   ``ray-circle``, 22
   ``singular-eigenfunction-eq46``, 21
   ``dd-curvilinear-scalar``, 20
   ``en-kernel-derivative``, 20
   ``kin-kernel-derivative``, 20
   ``peierls-rank-n-stability``, 20
   ``dd-cartesian-1d``, 19
   ``energy-condensation-rate-preservation``, 17
   ``direction-sampling``, 16
   ``discrete-measure-integrate``, 16
   ``fission-weight``, 16
   ``keff-cycle``, 16
   ``roulette-conservation``, 16
   ``roulette-prob``, 16
   ``transport-cartesian-2d``, 16
   ``peierls-greens-hollow-sph-architecture``, 15
   ``second-diff-general``, 15
   ``absorption-xs``, 14
   ``fission-source``, 14
   ``fixed-source-solve``, 14
   ``keff-update``, 14
   ``peierls-greens-annulus-architecture``, 14
   ``two-group-charpoly``, 14
   ``two-group-roots``, 14
   ``peierls-greens-slab-asym-architecture``, 13
   ``complementarity``, 12
   ``kinf-1g``, 12
   ``kinf-mg``, 12
   ``ld-cartesian-2d``, 12
   ``periodic-bc``, 12
   ``reciprocity``, 12
   ``ws-pitch``, 12
   ``loss-rep-scanmarch``, 11
   ``loss-rep-scanmarch-apply``, 11
   ``loss-rep-scanmarch-solve``, 11
   ``kll-1974-slab-flux``, 10
   ``kll-1974-sphere-flux``, 10
   ``peierls-greens-cylinder-architecture``, 10
   ``peierls-greens-slab-architecture``, 10
   ``peierls-vacuum-bc-flux``, 10
   ``peierls-vacuum-bc-row-sum-gate``, 10
   ``peierls-white-bc``, 10
   ``singular-eigenfunction-eq42``, 10
   ``dd-recurrence``, 9
   ``sn-curvilinear-trajectory-resolvent-crosscheck``, 9
   ``tau-m``, 9
   ``tau-p``, 9
   ``cp-inner-integral-antiderivative``, 8
   ``diffusion-coefficient``, 8
   ``energy-condensation-scattering-collapse``, 8
   ``bare-slab-buckling``, 7
   ``bare-slab-critical-equation``, 7
   ``bare-slab-eigenfunction``, 7
   ``dd-solve``, 7
   ``diffusion-M-matrix``, 7
   ``diffusion-back-substitution``, 7
   ``diffusion-exponential-branch``, 7
   ``diffusion-interface-matching``, 7
   ``diffusion-matching-matrix``, 7
   ``diffusion-mode-decomposition``, 7
   ``diffusion-operator``, 7
   ``diffusion-region-ode``, 7
   ``diffusion-spurious-root-validation``, 7
   ``diffusion-transcendental``, 7
   ``diffusion-trigonometric-branch``, 7
   ``gauss-legendre-visibility-cone``, 7
   ``inverse-as-operator``, 7
   ``bar-psi``, 6
   ``boyd-eq-45``, 6
   ``characteristic-ode``, 6
   ``kin-kernel-special-values``, 6
   ``ld-cartesian-1d``, 6
   ``peierls-greens-slab-T``, 6
   ``singular-eigenfunction-eq54``, 6
   ``sn-space-angle-separability``, 6
   ``cp-outer-integral-antiderivative``, 5
   ``dd-cartesian-2d``, 5
   ``en-kernel-special-values``, 5
   ``peierls-cyl-Gbc-3d-final``, 5
   ``peierls-greens-cylinder-T``, 5
   ``peierls-greens-cylinder-mr-homogeneous-reduction``, 5
   ``real-sh-discrete-orthogonality``, 5
   ``singular-eigenfunction-eq40``, 5
   ``streaming-equilibrium``, 5
   ``xs-interp``, 5
   ``dd-slab-scalar``, 4
   ``en-kernel-integral``, 4
   ``hebert-3-323``, 4
   ``ld-slab``, 4
   ``loss-rep-LpC``, 4
   ``macro-sum``, 4
   ``peierls-greens-V-alpha-1``, 4
   ``peierls-greens-V-alpha-2``, 4
   ``peierls-vacuum-bc-slab``, 4
   ``peierls-white-bc-slab``, 4
   ``phase-f-carlson-seed-source-driven``, 4
   ``phase-f-q-bar-twin-forms``, 4
   ``sn-mms-hetero-psi``, 4
   ``sn-mms-hetero-qext``, 4
   ``cp-escape-from-p-cell``, 3
   ``dd-mm-closure-constants``, 3
   ``delta-psi``, 3
   ``hetero-tolerance``, 3
   ``isotropic-source``, 3
   ``moc-keff-update``, 3
   ``moc-mms-psi-ref``, 3
   ``moc-mms-qext``, 3
   ``moc-wigner-seitz``, 3
   ``number-density``, 3
   ``peierls-greens-cylinder-mr-quadrature-convergence``, 3
   ``peierls-mg-operator``, 3
   ``peierls-vacuum-bc-cylinder``, 3
   ``peierls-vacuum-bc-sphere``, 3
   ``sigma-zero``, 3
   ``sn-homogenization-rate-preservation``, 3
   ``sn-mms-nonvacuum-psi``, 3
   ``sn-mms-nonvacuum-qext``, 3
   ``sn-space-angle-cross-term``, 3
   ``cp-second-difference-operator``, 2
   ``diffusion-mms``, 2
   ``harmonic-moment-projection``, 2
   ``hilbert-adjoint-equals-metric-times-S0``, 2
   ``peierls-greens-cylinder-mr-kinf``, 2
   ``peierls-greens-cylinder-mr-piecewise-tau``, 2
   ``peierls-greens-slab-asym-method-of-images``, 2
   ``peierls-slab-Gbc-mode``, 2
   ``peierls-slab-Pesc-mode``, 2
   ``roulette-restore``, 2
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
   ``sn-mms-cylindrical-psi``, 2
   ``sn-mms-cylindrical-qext``, 2
   ``sn-mms-nonvacuum-sph-psi``, 2
   ``sn-mms-nonvacuum-sph-qext``, 2
   ``sn-mms-p1-qext``, 2
   ``sn-mms-spherical-aniso-psi``, 2
   ``sn-mms-spherical-aniso-qext``, 2
   ``sn-mms-spherical-psi``, 2
   ``sn-mms-spherical-qext``, 2
   ``addition-theorem``, 1
   ``branching``, 1
   ``collision-estimator``, 1
   ``e1-decomposition``, 1
   ``majorant``, 1
   ``moment-projection-transpose-T``, 1
   ``nm1980-eq15-critical-condition``, 1
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
   ``peierls-greens-cylinder-mr-bounce-sum-piecewise``, 1
   ``peierls-greens-cylinder-mr-interface-continuity``, 1
   ``peierls-greens-cylinder-mr-trajectory-segments``, 1
   ``peierls-greens-cylinder-mr-wm72-vacuum``, 1
   ``peierls-greens-cylinder-trajectory``, 1
   ``peierls-greens-function-architecture``, 1
   ``peierls-greens-hollow-sph-impact-parameter-partition``, 1
   ``peierls-greens-hollow-sph-through-rank2``, 1
   ``peierls-greens-slab-asym-closure``, 1
   ``peierls-greens-slab-asym-resolvent``, 1
   ``peierls-greens-slab-trajectory``, 1
   ``peierls-greens-surface-fixed-point``, 1
   ``pi-r-equals-4pi-i``, 1
   ``real-spherical-harmonics``, 1
   ``sh-addition-theorem-reconstruction``, 1
   ``sh-space-metric``, 1
   ``si-spectral-rate``, 1
   ``sigT-computed``, 1
   ``singular-eigenfunction-eq5``, 1
   ``sn-mms-2d-2g-qext``, 1
   ``sn-mms-2d-psi``, 1
   ``sn-mms-2d-qext``, 1
   ``sn-mms-cylindrical-aniso-spatial-convergence``, 1
   ``sn-mms-p1-psi``, 1
   ``sn-mms-psi``, 1
   ``sn-mms-qext``, 1
   ``sn-mms-spherical-aniso-spatial-convergence``, 1
   ``splitting``, 1

Orphan equations
----------------

Equations with zero tests carrying ``@pytest.mark.verifies("label")``, excluding labels explicitly marked ``:vv-status: documented``. **112** of the testable equations found on theory pages are orphan.

- ``affine-bc-form``
- ``angular-windowing-moment-projection``
- ``apply-solve-cell-resolvent``
- ``apply-solve-denominator-inequality``
- ``apply-solve-neumann-expansion``
- ``apply-solve-neumann-series``
- ``apply-solve-parallel-identity``
- ``apply-solve-source-iteration-series``
- ``apply-solve-within-group-balance``
- ``bc-rank-n-tensor-decomposition``
- ``bc-tensor-decomposition``
- ``billiard-rank2-S``
- ``billiard-rank2-T``
- ``billiard-reflection-law``
- ``billiard-resolvent-neumann``
- ``billiard-transfer-operator``
- ``billiard-variant-alpha-rank1``
- ``dd-mm-angular-recurrence``
- ``dd-mm-scan-split``
- ``eigen-alpha-derivation``
- ``eigen-k-posing``
- ``eigen-resolvent``
- ``emission-spectrum-chi-mix``
- ``emission-spectrum-fission-source``
- ``emission-spectrum-simplex``
- ``fn-method-moment-space-AB-defs``
- ``fn-method-moment-space-bc-vacuum``
- ``fn-method-moment-space-fn-ansatz``
- ``fn-method-moment-space-galerkin-system``
- ``g-adjoint-block-metric``
- ``g-adjoint-derivation``
- ``g-adjoint-reciprocity``
- ``g-adjoint-sum-conjugation``
- ``g-adjoint-wrapper-action``
- ``galerkin-spectral-bte``
- ``galerkin-spectral-carlvik-integral``
- ``galerkin-spectral-eq4``
- ``galerkin-spectral-orthogonality``
- ``hebert-3-432-source``
- ``hebert-3-434``
- ``hebert-3-435``
- ``inflow-mask-discrete``
- ``ld-cartesian-2d-bilinear-coeffs``
- ``ld-cartesian-2d-face-bilinear-coeffs``
- ``ld-cartesian-2d-face-projection-coeff``
- ``ld-cartesian-2d-face-slot-shape``
- ``ld-cartesian-2d-projection-coeff``
- ``ld-ubld-cell-system``
- ``ld-ubld-d1-reduction``
- ``ld-ubld-divv-scale-free-kernel``
- ``ld-ubld-exact-on-bilinear``
- ``ld-ubld-kronecker-assembly``
- ``ld-ubld-kronecker-factors``
- ``ld-ubld-mass-weights``
- ``ld-ubld-moment-scan-source``
- ``ld-ubld-n-spatial-moments``
- ``ld-ubld-octant-moment-frame-signs``
- ``ld-ubld-pure-z-collision``
- ``ld-ubld-rule-of-three-collapse``
- ``ld-ubld-s2-s3-operators``
- ``ld-ubld-scale-free-invariants``
- ``ld-ubld-scattering-moment-lift``
- ``ld-ubld-slope-angular-reduction``
- ``ld-ubld-unified-moment-residual``
- ``ld-ubld-weak-form``
- ``loss-rep-affine``
- ``loss-rep-affine-cell``
- ``loss-rep-facewise-separable``
- ``loss-rep-leaf-sum``
- ``loss-rep-removal-sigma``
- ``loss-rep-scanmarch-apply-residual``
- ``loss-rep-scanmarch-solve-affine``
- ``ordinate-partition-inflow-outflow``
- ``phase-c-cell-update``
- ``phase-c-streaming-spherical``
- ``phase-c-wdd-oscillation``
- ``phase-c-wdd-recurrence``
- ``phase-f-q-1d-decomposition``
- ``phase-f-source-eq-sigt-phi0``
- ``pole-mm-recurrence``
- ``richardson-diffusion``
- ``scattering-spectral-theorem``
- ``si-jacobi-fixed-point``
- ``si-sigma-r-fold-mismatch``
- ``si-within-group-operator-eq``
- ``sn-affine-outgoing-face-reconstruction-eq``
- ``sn-axis-widths``
- ``sn-cell-flatten-roundtrip``
- ``sn-p1-cylinder-hand-ref``
- ``sn-p1-sphere-hand-ref``
- ``sn-pole-cell-shell-average``
- ``sn-streaming-reciprocity``
- ``sn-tau-mm-raw``
- ``sn-within-group-system``
- ``spatial-moment-append-policy``
- ``spatial-moment-kronecker-order``
- ``spatial-moment-space-size``
- ``spectrum-case-eigenfunction-equation``
- ``spectrum-case-eigenfunction-explicit``
- ``spectrum-continuum-eigenfunction``
- ``spectrum-dispersion-relation``
- ``spectrum-expansion-theorem``
- ``spectrum-full-decomposition``
- ``spectrum-transport-equation``
- ``spectrum-x-function``
- ``streaming-action-cell-balance``
- ``trace-sign-predicate``
- ``two-moment-angular``
- ``two-moment-carrier-space``
- ``two-moment-spatial``
- ``two-moment-tensor-product``
- ``vacuum-legacy-vs-trace-correct``

Documented-only equations
-------------------------

Theory labels marked ``.. vv-status: <label> documented`` in their RST source. These are excluded from the orphan-equation gate because they are either definitional (no single implementing function — e.g. ``boltzmann``), describe a module whose Python port does not yet exist (e.g. the thermal-hydraulics / fuel-behaviour / reactor-kinetics equations), or have a deliberately deferred test paired with a tracking issue. **280** labels carry the directive. See ``docs/testing/architecture.rst``:ref:`vv-status-documented` for the full taxonomy.

- ``affine-contraction-ratio``
- ``affine-torsor-algebra``
- ``affine-true-error``
- ``affine-typed-residual-eq``
- ``angular-windowing-aniso-factoring``
- ``angular-windowing-moment-iterate``
- ``apply-distributes``
- ``bailey-dome-recursion``
- ``bc-extraction-block-matrix``
- ``bc-extraction-direct-sum-state``
- ``bc-extraction-loss-operator``
- ``bc-extraction-two-hat-closed-sums``
- ``bc-extraction-two-residuals``
- ``bc-extraction-variadic-matvec``
- ``bc-extraction-within-group-decomposition``
- ``bc-rank-n-as-sum-of-products``
- ``bessel-wronskian``
- ``bickley-integral``
- ``boltzmann``
- ``bundle-measure-disintegration``
- ``burst-criterion``
- ``c-in-remapping``
- ``carrier-grid-cell``
- ``carrier-grid-interchange-witness``
- ``carrier-grid-operator-typing``
- ``case-dispersion-function``
- ``clad-heat``
- ``conservative-form``
- ``coolant-energy``
- ``coolant-feedback``
- ``coolant-rate``
- ``creep-rate``
- ``dd-2d-balance-form``
- ``diagonal-operator-action``
- ``diffusion-boundary-closure``
- ``diffusion-interior-conductance``
- ``diffusion-operator-family``
- ``diffusion-partial-current-dictionary``
- ``discrete-measure-definition``
- ``discrete-measure-partition``
- ``discrete-measure-pushforward``
- ``doppler-feedback``
- ``e1-small-tau-expansion``
- ``eigen-standard-form``
- ``energy-condensation-balance``
- ``energy-condensation-chi-collapse``
- ``energy-condensation-coarse-flux``
- ``energy-condensation-counting-measure``
- ``energy-condensation-fine-rate``
- ``energy-condensation-lethargy-overlap``
- ``energy-condensation-matrix-collapse``
- ``energy-condensation-nested-subset``
- ``energy-condensation-overlap-fraction``
- ``energy-condensation-partition-of-unity``
- ``energy-condensation-representative-spectrum``
- ``energy-condensation-vector-collapse``
- ``fb-bc4-displacement``
- ``fb-clad-strain``
- ``fb-fuel-heat``
- ``fb-fuel-strain``
- ``fb-swelling``
- ``fission-as-dyad``
- ``fn-Fk-closed-forms``
- ``fn-Fk-integration-by-parts``
- ``fn-product-simpson-weights``
- ``fn-slab-B-long-division``
- ``fn-slab-collocation``
- ``fn-unified-matrix-entry``
- ``fn-x-function``
- ``fuel-heat``
- ``fuel-rate``
- ``funk-hecke-eigenvalue``
- ``g-adjoint-definition``
- ``galerkin-construction``
- ``galerkin-pair``
- ``galerkin-self-adjoint``
- ``gap-closure-event``
- ``gap-conductance``
- ``gas-pressure``
- ``green-neumann-series``
- ``green-splitting-iteration``
- ``green-true-residual``
- ``group-flux``
- ``group-xs``
- ``hilbert-adjoint-equals-metric-times-S0``
- ``in-scatter-full-contraction``
- ``integral-kernel-category``
- ``keff-as-integrated-rates``
- ``kll-1974-slab-phi``
- ``kll-1974-sphere-phi``
- ``matrix-functor-out``
- ``matrix-inverse-direct-residual``
- ``matrix-inverse-materialise``
- ``maxwellian``
- ``mm-half-grid-recurrence``
- ``mode-conservation-target``
- ``moment-projection-transpose-T``
- ``morel-montry-clamp``
- ``multiplication-operator-action``
- ``multiplication-operator-embedding``
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
- ``peierls-greens-cylinder-mr-bounce-sum-piecewise``
- ``peierls-greens-cylinder-mr-piecewise-tau``
- ``peierls-greens-cylinder-mr-trajectory-segments``
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
- ``peierls-mg-fission-source-local``
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
- ``product-solve-reroute``
- ``production-rate-functional``
- ``quadrature-selection-criterion``
- ``reaction-rate-kinf-oracle``
- ``real-sh-addition-theorem``
- ``real-sh-discrete-orthogonality``
- ``real-sh-l0``
- ``real-sh-l1``
- ``real-sh-l2plus``
- ``resolvent-similarity``
- ``scattering-as-tensor-product-sum``
- ``scattering-carrier-grid``
- ``scattering-zonal-kernel``
- ``sh-addition-theorem-reconstruction``
- ``sh-space-metric``
- ``si-within-group-fixed-point``
- ``sigs-convention``
- ``singular-eigenfunction-eq5``
- ``sn-coupled-pole-mu-level-invariant-eq``
- ``sn-err-058-coupled-pole-continuity``
- ``sn-err-058-edge-extrapolation``
- ``sn-err-058-proxy-source``
- ``sn-homogenization-balance``
- ``sn-homogenization-bilinear``
- ``sn-homogenization-chi-collapse``
- ``sn-homogenization-coarse-space``
- ``sn-homogenization-fine-rate``
- ``sn-homogenization-frame-projector``
- ``sn-homogenization-matrix-collapse``
- ``sn-homogenization-metric-fold``
- ``sn-homogenization-normal-equations``
- ``sn-homogenization-production-weight``
- ``sn-homogenization-radon-nikodym``
- ``sn-homogenization-region-flux``
- ``sn-homogenization-scatter-rate``
- ``sn-homogenization-test-functions``
- ``sn-homogenization-vector-collapse``
- ``sn-keff-update``
- ``sn-mms-nonvacuum-psi``
- ``sn-mms-nonvacuum-sph-psi``
- ``sn-space-angle-cross-term``
- ``sn-space-angle-separability``
- ``solve-does-not-distribute``
- ``sood-eq18-1g-balance``
- ``sood-eq19-kinf-1g``
- ``sood-eq20-kinf-1g-c-form``
- ``sood-eq25-2g-matrix``
- ``sood-eq29-kinf-2g-no-upscatter``
- ``sood-eq32-phi-ratio``
- ``sood-eq76-kinf-mg``
- ``streaming-action-pure-l``
- ``streaming-as-tensor-product-sum``
- ``streaming-inverse-direct-sum``
- ``streaming-pn-recurrence``
- ``subgroup-of-o3-containment``
- ``sum-of-tensor-products``
- ``tensor-product-action``
- ``tensor-product-adjoint-distributivity``
- ``tensor-product-axis-wise-composition``
- ``tensor-product-inverse``
- ``trace-half-decomposition``
- ``transport-equation``
- ``vacuum-bc``
- ``wave-t-cell-balance-three-terms``
- ``wave-t-ma-q1-master-condition``
- ``wave-t-mspat-curvilinear-subtraction``
- ``wavefront-cochain-biproduct``
- ``wavefront-cochain-biproduct-laws``
- ``wavefront-cochain-primal``
- ``wdd-forward-recurrence``
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
   ``ERR-003``, 2
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
   ``ERR-025``, 5
   ``ERR-026``, 98
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
   ``ERR-039``, 8
   ``ERR-040``, **0 (MISSING)**
   ``ERR-041``, **0 (MISSING)**
   ``ERR-042``, **0 (MISSING)**
   ``ERR-043``, 8
   ``ERR-044``, 4
   ``ERR-045``, **0 (MISSING)**
   ``ERR-046``, 8
   ``ERR-047``, **0 (MISSING)**
   ``ERR-048``, 51
   ``ERR-049``, 15
   ``ERR-050``, 4
   ``ERR-051``, **0 (MISSING)**
   ``ERR-052``, 1
   ``ERR-053``, 7
   ``ERR-054``, 2
   ``ERR-055``, **0 (MISSING)**
   ``ERR-056``, 7
   ``ERR-057``, 1
   ``ERR-058``, 5
   ``ERR-059``, 4
   ``ERR-060``, 2
   ``ERR-061``, 3
   ``ERR-062``, 1
   ``ERR-063``, 3
   ``ERR-064``, 2
   ``ERR-065``, 1

Unmarked tests
--------------

**10 tests** have no V&V level marker.
This is a gap — every test in the tree should carry either
a physics-ladder marker (``l0``..``l3``) or the orthogonal
``foundation`` marker (``@pytest.mark.foundation``) for
tests that verify software invariants rather than physics
equations. See ``docs/testing/architecture.rst``
:ref:`vv-foundation-tests` for the taxonomy.

.. csv-table::
   :header: File, Unmarked tests
   :widths: 60, 10

   ``tests/derivations/test_atkinson_product_nystrom.py``, 3
   ``tests/derivations/test_fn_projection_vs_kll_flux.py``, 2
   ``tests/sn/sweep/core/test_cell_kernel_batch.py``, 2
   ``tests/derivations/test_fn_la13511_slab_reflected.py``, 1
   ``tests/sn/spatial/test_linear_discontinuous.py``, 1
   ``tests/sn/verification/mms/test_mms_ld_2d.py``, 1

