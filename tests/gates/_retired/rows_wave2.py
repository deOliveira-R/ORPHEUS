"""The rows of wave 2 (WF2.4, the typed-law vocabulary and its three homes), ported from ``tests/gates/test_wave2_retirement.py``.

60 ``RetiredName`` rows (the retiring names with no surviving homonym; each reads 0 identifier hits in the five trees, the measured
count of the port, 60 of 60 green), 2 ``Retired`` rows (the two properties ``geometry_map`` and ``response_kernel`` of the base class
``BoundaryTraceLaw``, bound by definition site because the bare words are generic locals, decision D5) and 3 ``RetiredModule`` rows
(the three retired homes; each also asserts that its file and package directory are absent, decision E-D6). ``KEPT = ()``: the
survivors that the old file's ``_HOMONYMS`` exempted exempt nothing under the narrowed reading (decision R2).
"""

from __future__ import annotations

from tests.gates._retired import Kept, Retired, RetiredModule, RetiredName, Row

ROWS: tuple[Row, ...] = (
    # the 60 bare names: each reads 0 identifier hits in the five trees at the port (60 of 60 green)
    RetiredName('AlbedoBoundary'),
    RetiredName('BC'),
    RetiredName('BC_REGISTRY'),
    RetiredName('BOUNDARY_OPERATOR_REGISTRY'),
    RetiredName('BoundaryGeometryMap'),
    RetiredName('BoundaryRealizer'),
    RetiredName('BoundaryResponseKernel'),
    RetiredName('BoundaryTraceLaw'),
    RetiredName('ConstantInflowSource'),
    RetiredName('DiffusionBoundaryRealizer'),
    RetiredName('FaceLaws'),
    RetiredName('InflowSourceSpec'),
    RetiredName('IsotropicReturn'),
    RetiredName('LambertianReemission'),
    RetiredName('LawNode'),
    RetiredName('LawScaled'),
    RetiredName('LawSum'),
    RetiredName('NoSource'),
    RetiredName('PairedDeck'),
    RetiredName('PeriodicBoundary'),
    RetiredName('PrescribedInflow'),
    RetiredName('ReemissionClosure'),
    RetiredName('ReflectiveBoundary'),
    RetiredName('SNBoundaryRealizer'),
    RetiredName('ScalarResponse'),
    RetiredName('SelfPairedDeck'),
    RetiredName('SpecularReemission'),
    RetiredName('SpecularReturn'),
    RetiredName('TAG_REGISTRY'),
    RetiredName('VacuumAppliedToOutgoingTraceError'),
    RetiredName('VacuumInflow'),
    RetiredName('WhiteBoundary'),
    RetiredName('ZeroFluxBoundary'),
    RetiredName('_BoundBoundaryOperator'),
    RetiredName('_admit_as_operand'),
    RetiredName('_axis_cosines'),
    RetiredName('_checked_angular_average'),
    RetiredName('_deck_kernel'),
    RetiredName('_law_from_tag'),
    RetiredName('_law_of'),
    RetiredName('_mirror_motion'),
    RetiredName('_narrowed_zero_operator'),
    RetiredName('_specular_crossed_diagnosis'),
    RetiredName('_specular_pairing'),
    RetiredName('_wall_of'),
    RetiredName('assert_realizable'),
    RetiredName('assert_specular_pairing_involutive'),
    RetiredName('assert_specular_pairing_maps_inflow_to_outflow'),
    RetiredName('assert_specular_pairing_measure_preserving'),
    RetiredName('assert_specular_pairing_valid'),
    RetiredName('inflow_indices'),
    RetiredName('law_permutes_ordinates'),
    RetiredName('outflow_indices'),
    RetiredName('parse_boundary_law'),
    RetiredName('parse_boundary_laws'),
    RetiredName('realize_boundary_law'),
    RetiredName('realize_recursively'),
    RetiredName('resolve_boundary_conditions'),
    RetiredName('specular_albedo'),
    RetiredName('to_alpha'),
    # the two properties of the retired base class: 0 definitions, 0 setter or import uses (2 of 2 green)
    Retired('orpheus/geometry/boundary/_base.py', 'BoundaryTraceLaw.geometry_map'),
    Retired('orpheus/geometry/boundary/_base.py', 'BoundaryTraceLaw.response_kernel'),
    # the three retired homes: 0 imports and no file or package directory (3 of 3 green)
    RetiredModule('orpheus.geometry.boundary'),
    RetiredModule('orpheus.transport.boundary_realization'),
    RetiredModule('orpheus.sn.boundary.bound_compat'),
)
KEPT: tuple[Kept, ...] = ()
