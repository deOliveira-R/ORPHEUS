"""The rows of sub-unit dead-A of unit ``dead`` (specification ``ontology_w0_spec.md`` Appendix D-A.7): measure, quadrature selector,
geometry and descent.

``ROWS``: 13 ``Retired`` rows, plus ``QuadratureSpec.structural_flags``, which is no row of ``reviewed.json`` but the ruled closure
member of row 73 (3.0 item 3), and 11 ``Moved`` rows, each with its concrete home and the name it carries there. ``KEPT``: the 30
rows the unit does not touch (the 26 geometry and measure rows, and rows 134 to 137, the four inverters kept as catalogue data for
W8a, 3.0 item 2), asserted still defined, so that a builder cannot retire a kept symbol by accident.

Every row is defined on the base tree, so each ``Retired`` and ``Moved`` row is red until the unit lands; ``KEPT`` is green on both.
"""

from __future__ import annotations

from tests.gates._retired import Kept, Moved, Retired, Row

ROWS: tuple[Row, ...] = (
    Moved('orpheus/geometry/chart.py', 'DirectionDomain.tangencies', 'tests/gates/geometry/_geometry_oracles.py', 'tangencies'),  # row 6; free function tangencies(direction_domain, level)
    Retired('orpheus/geometry/chord.py', 'Chord.slot_start'),
    Retired('orpheus/geometry/coord.py', 'compute_volumes_1d'),
    Moved('orpheus/geometry/manifold.py', 'UNIT_INTERVAL', 'tests/gates/geometry/_manifold_fixtures.py', 'UNIT_INTERVAL'),  # row 27; a value
    Retired('orpheus/geometry/manifold.py', 'barycentre'),
    Retired('orpheus/geometry/structured_geometry.py', 'StructuredGeometry.domain_extent_cm'),
    Retired('orpheus/geometry/symmetry.py', 'SubgroupOfO3.identity_component'),
    Moved('orpheus/geometry/symmetry.py', '_group_elements', 'tests/gates/geometry/_geometry_oracles.py', '_group_elements'),  # row 32; free function _group_elements(tag)
    Moved('orpheus/geometry/symmetry.py', 'OrbitCertificate.stabilizer_order', 'tests/gates/geometry/_geometry_oracles.py', 'stabilizer_order'),  # row 33; free function stabilizer_order(orbit_certificate)  # was a property
    Moved('orpheus/geometry/transformation.py', 'RigidMotion.element_order', 'tests/gates/geometry/test_transformation.py', 'element_order'),  # row 34; free function element_order(rigid_motion, *, atol=1e-9)
    Moved('orpheus/geometry/transformation.py', 'RigidMotion.permutes', 'tests/gates/geometry/_geometry_oracles.py', 'permutes'),  # row 36; free function permutes(rigid_motion, points, *, atol)
    Retired('orpheus/numerics/basis/descent.py', 'Descent.discriminator'),
    Retired('orpheus/numerics/basis/descent.py', 'Descent.descending_real_slots'),
    Retired('orpheus/numerics/basis/descent.py', 'Descent.is_isomorphism'),
    Retired('orpheus/numerics/elementary_operators.py', 'DiagonalOperator.from_measure'),
    Moved('orpheus/numerics/generating_measure.py', 'HERMITE', 'tests/gates/numerics/_measure_fixtures.py', 'HERMITE'),  # row 54; a value
    Moved('orpheus/numerics/generating_measure.py', 'laguerre', 'tests/gates/numerics/_measure_fixtures.py', 'laguerre'),  # row 55; free function laguerre(a=0.0)
    Moved('orpheus/numerics/index_set.py', 'ENERGY', 'tests/gates/numerics/_measure_fixtures.py', 'ENERGY'),  # row 56; a value
    Moved('orpheus/numerics/measure.py', 'equispaced', 'tests/gates/numerics/_measure_fixtures.py', 'equispaced'),  # row 60; free function equispaced(a, b, n)
    Moved('orpheus/numerics/mesh_free_function.py', 'Symbolic.from_srepr', 'tests/gates/numerics/_measure_fixtures.py', 'from_srepr'),  # row 61; free function from_srepr(srepr, sympy_version=None), its body names `Symbolic`  # was a classmethod: the class parameter is dropped (decision D2)
    Retired('orpheus/numerics/quadrature/registry.py', 'SelectionLog.summary'),
    Retired('orpheus/numerics/quadrature/registry.py', 'select_quadrature'),
    Retired('orpheus/numerics/quadrature/registry.py', 'SelectionLog'),
    Retired('orpheus/numerics/quadrature/registry.py', 'QuadratureSelectionError'),
    Retired('orpheus/numerics/quadrature/registry.py', 'QuadratureSpec.structural_flags'),  # not a row: the ruled closure member of row 73 (3.0 item 3, 8.24)
)

KEPT: tuple[Kept, ...] = (  # Kept(file, qualname, keeper): the keeper is the wave or the consumer that keeps the row; asserted still defined
    Kept('orpheus/geometry/chart.py', 'Chart.measure_density', 'derivations/ (A2): characteristic/basis.py:285, closure.py:284'),  # row 3
    Kept('orpheus/geometry/chart.py', 'Chart.directions_at', 'derivations/ (A2): characteristic/reading.py:176,202,276'),  # row 4
    Kept('orpheus/geometry/chart.py', 'Chart.line_domain', 'derivations/ (A2): characteristic/assembly.py, lines.py, reading.py'),  # row 5
    Kept('orpheus/geometry/chart.py', 'RadialImage.orbit_coordinate_at', 'derivations/ (A2) via Chord/Rays.orbit_coordinate_at'),  # row 7
    Kept('orpheus/geometry/chart.py', 'AxialImage.orbit_coordinate_at', 'as row 7'),  # row 8
    Kept('orpheus/geometry/chord.py', 'ConcentricPartition.region_at', "W6, concept row 'Point location, chord' (#596: ConcentricPartition is the one locator)"),  # row 9
    Kept('orpheus/geometry/chord.py', 'ConcentricPartition.refined_at', 'derivations/ (A2): characteristic/reading.py:204,280'),  # row 10
    Kept('orpheus/geometry/chord.py', 'ConcentricPartition.coarsening', 'derivations/ (A2): characteristic/basis.py:143, reading.py:206, transport.py:679'),  # row 11
    Kept('orpheus/geometry/chord.py', 'ConcentricPartition.rays_from', 'derivations/ (A2): peierls geometry.py, fn sphere, characteristic'),  # row 12
    Kept('orpheus/geometry/chord.py', 'Chord.orbit_coordinate_at', 'derivations/ (A2)'),  # row 14
    Kept('orpheus/geometry/chord.py', 'Chord.rays_at', 'derivations/ (A2): characteristic/transport.py:453'),  # row 15
    Kept('orpheus/geometry/chord.py', 'Chord.rays_at_closest_approach', 'derivations/ (A2): characteristic/transport.py:456'),  # row 16
    Kept('orpheus/geometry/chord.py', 'Chord.line_integral', 'derivations/ (A2): 6 reads'),  # row 17
    Kept('orpheus/geometry/chord.py', 'Chord.cosine_at', 'derivations/ (A2): peierls geometry.py (6)'),  # row 18
    Kept('orpheus/geometry/chord.py', 'Rays.lengths_to', 'derivations/ (A2): peierls geometry.py:910,1181,1263'),  # row 19
    Kept('orpheus/geometry/chord.py', 'Rays.to_wall', 'derivations/ (A2): 6 reads'),  # row 20
    Kept('orpheus/geometry/chord.py', 'Rays.orbit_coordinate_at', 'derivations/ (A2)'),  # row 21
    Kept('orpheus/geometry/chord.py', 'Rays.advanced', 'derivations/ (A2): characteristic/transport.py:459'),  # row 22
    Kept('orpheus/geometry/coord.py', 'CoordSystem.measure_density', 'derivations/ (A2) through Chart.measure_density'),  # row 23
    Kept('orpheus/geometry/deck.py', 'ISOMETRY_CLASSES_BY_DIMENSION', 'W6 (declared deferral table, read with DEFERRED_DECKS)'),  # row 25
    Kept('orpheus/geometry/deck.py', 'DEFERRED_DECKS', "W6: SCOPE-BOUNDARY table named by production's DeckError message (geometry/deck.py:293)"),  # row 26
    Kept('orpheus/geometry/structured_geometry.py', 'StructuredGeometry.is_hollow', 'derivations/ (A2): common/reference_body.py:164,188'),  # row 29
    Kept('orpheus/geometry/transformation.py', 'RigidMotion.seated_at', "W6, concept row 'Symmetry group, translations included' (a subgroup of E(d))"),  # row 35
    Kept('orpheus/numerics/measure.py', 'DiscreteMeasure.singular_set_under', "W8a, concept row 'Discrete measure' (the invariance half of the weld)"),  # row 57
    Kept('orpheus/numerics/measure.py', 'DiscreteMeasure.symmetry_groups', 'W8a, as row 57'),  # row 58
    Kept('orpheus/numerics/measure.py', 'BundleMeasure', "W8a, concept row 'Discrete measure' ('widened to differing supports')"),  # row 59
    Kept('orpheus/numerics/quadrature/registry.py', '_gl1d_invert', 'W8a: a value of the required field degree_of_exactness_for of the kept catalogue (3.0 item 2)'),  # row 134
    Kept('orpheus/numerics/quadrature/registry.py', '_lebedev_invert', 'as row 134'),  # row 135
    Kept('orpheus/numerics/quadrature/registry.py', '_ls_sn_invert', 'as row 134'),  # row 136
    Kept('orpheus/numerics/quadrature/registry.py', '_product_invert', 'as row 134'),  # row 137
)
