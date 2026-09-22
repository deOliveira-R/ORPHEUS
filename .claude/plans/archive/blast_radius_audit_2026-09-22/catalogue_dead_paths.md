# Error-catalogue test paths with NO successor by basename (measured 2026-09-22)

Predicate: `tests/…\.py` paths cited in `docs/theory/verification/error_catalog.rst` that do not exist and whose basename matches no file under `tests/` (`pathlib.rglob`). The 16 one-successor moves were re-pointed at `cc51e027`. Each row: the dead path, then the ERR entry and line of every citation.

| dead path | cited by (entry:line) |
|---|---|
| `tests/derivations/test_peierls_slab_multiregion.py` | ERR-033:2322 |
| `tests/numerics/test_projection_operators.py` | ERR-039:3204, ERR-039:3257, ERR-051:4570, ERR-051:4592 |
| `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` | ERR-026:1331, ERR-026:1460, ERR-026:1523, ERR-026:1524 |
| `tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py` | ERR-026:1456 |
| `tests/sn/operators/test_starting_direction_metric.py` | ERR-067:5224 |
| `tests/sn/spatial/test_boundary_face_flux.py` | ERR-026:1486 |
| `tests/sn/spatial/test_pole_angular_closure.py` | ERR-026:1446 |
| `tests/sn/test_cartesian.py` | ERR-025:1128 |
| `tests/sn/test_cylindrical.py` | ERR-055:4906, ERR-055:4931 |
| `tests/sn/test_invertible_operator.py` | ERR-049:4310, ERR-049:4320, ERR-049:4340, ERR-050:4523 |
| `tests/sn/test_phase_c_mms.py` | ERR-026:1218 |
| `tests/sn/test_snstreamingoperator.py` | ERR-026:1389, ERR-026:1629 |
| `tests/sn/test_spherical.py` | ERR-055:4906, ERR-055:4931 |
| `tests/sn/test_sweep_operator_inconsistency.py` | ERR-026:1212, ERR-026:1305, ERR-026:1441 |
