# Catalogue test-FUNCTION citations that do not resolve (AST, 2026-09-22)

Predicate: `tests/…py::A[::B]` citations in `error_catalog.rst` whose file exists and whose named def/class/module-level name does not (a trailing `{`, `*` or `[` is read as a prefix pattern). Positive control: ERR-087's catcher resolves.

**12 of 77** do not resolve:

- line 752: `tests/geometry/test_geometry.py::TestZoneSubdivision::test_equal_volume_single_zone`
- line 753: `tests/geometry/test_geometry.py::TestZoneSubdivision::test_equal_volume_multi_zone`
- line 790: `tests/geometry/test_geometry.py::TestZoneSubdivision::test_equal_volume_`
- line 839: `tests/moc/test_ray_tracing.py::test_degenerate_corner_ray`
- line 4938: `tests/sn/sweep/core/test_ordinate_scan.py::test_ordinate_scan_small_attenuation`
- line 4964: `tests/sn/sweep/core/test_ordinate_scan.py::test_pair_monoid_associativity`
- line 5182: `tests/transport/spatial/test_ld_ubld_symbolic.py::test_d2_exact_on_bilinear`
- line 5182: `tests/transport/spatial/test_ld_ubld_primitive.py::test_d2_exact_on_bilinear`
- line 5182: `tests/transport/spatial/test_linear_discontinuous.py::test_d2_assembled_matrices_match_symbolic`
- line 5280: `tests/moc/test_verification.py::test_n2n_1g_analytical_keff`
- line 6850: `tests/numerics/test_basis_domain.py::test_e1`
- line 6866: `tests/numerics/test_quadrature_directional.py::test_q8_4_the_1d_lift_is_still_a_FICTION_and_says_so`
