r"""Two owners, one space: each ``SNProblem`` phase-space mint equals the raw constructor applied to the problem's own
ingredients (layer DAG Phase 4, ruling 24).

Since Phase 4 every mint is ONE call to ``orpheus.transport.phase_spaces`` (or the scheme's ``trial_space``), so these
rows are DEMOTED with respect to the RULE (``retirement-audit`` D.16): both sides run the same function, and a
mutation of that function leaves the eight mint rows green (``[M]`` arm ``bulk_uniform_ordinate_weights``: 0 of 8; it
reds the 3 moment-entry rows, whose frame side reads the mutated ordinate axis); the rule's reference is
``tests/gates/transport/test_phase_spaces.py``. What they still test is the ROUTE: the problem hands the
constructor ITS quadrature and ITS scheme, and its material-mesh data equals a standalone ``MaterialMesh`` built from
the same mesh and materials. First reds: ``sn_trace_reads_a_fixed_scheme`` (the trace mint built with a Diamond
Difference scheme whatever the problem binds: red on the LD rectangle, where the face layout carries the moment tail)
and ``sn_bulk_reads_a_fixed_quadrature``.

The rows also pin the one identity the hub keeps that a raw constructor does not: a slopeless problem's trial space
IS its bulk space (one cached instance).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.sn.problem import SNProblem
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.phase_spaces import angular_bulk_space, angular_trace_space
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from tests.gates._exterior import conditions_for
from tests.gates.transport import _phase_space_cases as cases
from tests.gates.transport._materials import CURVED_REFLECTIVE

pytestmark = pytest.mark.foundation

_CURVED = ("sphere", "hollow_sphere", "cylinder", "hollow_cylinder")
_ROWS = [(g, s) for g in cases.GEOMETRIES for s in ("DD",) + (() if g in _CURVED else ("LD",))]


def _pose(geometry: str, scheme: str) -> tuple[SNProblem, MaterialMesh]:
    material_mesh = cases.GEOMETRIES[geometry]()
    mesh = cases.legacy_mesh(material_mesh)
    conditions = conditions_for(mesh, {0: WallCondition.vacuum()} if geometry not in _CURVED else CURVED_REFLECTIVE)
    problem = SNProblem(
        mesh, cases.QUADRATURES[geometry](), material_mesh.materials,
        scheme=LinearDiscontinuous() if scheme == "LD" else DiamondDifference(), boundary_conditions=conditions,
    )
    return problem, MaterialMesh(mesh, material_mesh.materials)


@pytest.mark.parametrize("geometry,scheme", _ROWS, ids=[f"{g}-{s}" for g, s in _ROWS])
@pytest.mark.rests_on(
    "tests/gates/transport/test_phase_spaces.py::test_the_angular_bulk_is_the_ordinate_axis_on_the_meshs_own_axes",
    "tests/gates/transport/test_phase_spaces.py::test_the_angular_trace_is_the_numerics_trace_over_the_hand_layout",
    "tests/gates/transport/spatial/test_trial_space.py::test_a_slopeless_closure_returns_the_space_itself",
    "tests/gates/transport/spatial/test_trial_space.py::test_a_multi_moment_closure_appends_the_moment_axis_with_its_mass",
)
def test_the_problems_mints_are_the_constructors_on_its_own_ingredients(geometry: str, scheme: str) -> None:
    problem, standalone = _pose(geometry, scheme)
    bulk = angular_bulk_space(standalone, problem.quad)
    trial = problem.scheme.trial_space(bulk, standalone.axes)
    trace = angular_trace_space(standalone, problem.quad, problem.scheme)
    for name, mint, built in (
        ("angular_bulk_space", problem.angular_bulk_space, bulk),
        ("angular_trial_space", problem.angular_trial_space, trial),
        ("angular_trace", problem.angular_trace, trace),
        ("full_field_space", problem.full_field_space, FullFieldSpace.from_blocks(trial, trace)),
    ):
        if mint != built or hash(mint) != hash(built):
            pytest.fail(f"{geometry}/{scheme}: SNProblem.{name} differs from the constructor on the problem's own ingredients")
    if scheme == "DD" and problem.angular_trial_space is not problem.angular_bulk_space:
        pytest.fail(f"{geometry}: a slopeless problem's trial space is not its bulk space (one cached instance)")


@pytest.mark.parametrize("geometry,scheme", [("slab", "DD"), ("cart2d", "LD")], ids=["slab-DD", "cart2d-LD"])
def test_the_carrier_keyed_packing_is_the_space_keyed_packing_on_the_problems_trace(geometry: str, scheme: str) -> None:
    problem, _ = _pose(geometry, scheme)
    rng = np.random.default_rng(3)
    arrays = {key: rng.standard_normal(slot.shape) for key, slot in problem.angular_trace.layout.faces.items()}
    by_carrier = AngularBoundaryFlux.from_face_arrays(problem, arrays)
    by_space = AngularBoundaryFlux.from_face_arrays_on(problem.angular_trace, arrays)
    if by_carrier.space is not problem.angular_trace or by_space.space is not problem.angular_trace:
        pytest.fail(f"{geometry}/{scheme}: a packed field does not hold the problem's trace object")
    np.testing.assert_array_equal(np.asarray(by_carrier.values), np.asarray(by_space.values))


@pytest.mark.parametrize("geometry,scheme", [("slab", "LD"), ("cart2d", "LD"), ("slab", "DD")], ids=["slab-LD", "cart2d-LD", "slab-DD"])
def test_the_carrier_keyed_moment_entries_are_the_space_keyed_spellings(geometry: str, scheme: str) -> None:
    """The two remaining carrier-keyed entries of census F1, against their space-keyed spellings on the problem's own
    ingredients: the composer ``BulkField.compose_spatial_moments(space, problem, width)`` is the scheme's
    ``trial_space(space, axes)`` at the scheme's own width (``space`` itself at width 1), and the factory
    ``HarmonicMomentFlux.zeros_for_problem_and_L(problem, L, spatial_moments=width)`` holds the frame's
    ``moment_space_on`` of the (widened) angular bulk."""
    from orpheus.transport.fields._bases import BulkField
    from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
    from orpheus.transport.frames import HarmonicFrame

    problem, standalone = _pose(geometry, scheme)
    width = problem.scheme.spatial_basis_per_axis
    bulk = angular_bulk_space(standalone, problem.quad)
    widened = problem.scheme.trial_space(bulk, standalone.axes)
    composed = BulkField.compose_spatial_moments(problem.angular_bulk_space, problem, width)
    if composed != widened:
        pytest.fail(f"{geometry}/{scheme}: compose_spatial_moments at width {width} is not the scheme's trial space")
    for L in (1, 2):
        field = HarmonicMomentFlux.zeros_for_problem_and_L(problem, L, spatial_moments=width)
        expected = HarmonicFrame.for_space(widened, L).moment_space_on(widened)
        if field.space != expected:
            pytest.fail(f"{geometry}/{scheme} L={L}: the factory's moment space is not the frame's on the widened bulk")
