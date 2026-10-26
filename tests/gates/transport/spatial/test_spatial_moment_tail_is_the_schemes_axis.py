r"""The moment product's spatial-moment TAIL is the scheme's own axis — one
spelling of the factor on every side (CS4c step 6 item 6.2c-iii, 2026-09-08).

Until this item the tree spelled the within-cell spatial-moment factor two
ways: the widened ANGULAR and SCALAR spaces carried the scheme-owned MODAL
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.moment_axis`
(mass-weighted — its weights are the cell basis's Gram; the field
layer's composer appends it by label), while the harmonic-MOMENT product
appended a separate Euclidean, axes-less ``SpatialMomentSpace`` class — so a
widened moment space was axes-less (its ``*`` product left the axis arm), its
tail carried NO mass, and the frame's derivation had to DROP the angular
space's tail axis and re-append the class to stay ``(name, shape)``-equal to
the hub's. The carrier's hub now composes its tail THROUGH the fields' own
composer and the frame keeps the angular space's own axes, so:

* the widened moment product is axis-built, carries the scheme's axis as its
  tail, and equals the frame's derivation structurally (ruling O-5 at width
  2 as at width 1);
* the moment field's norm on a widened iterate carries the cell mass on its
  tail exactly as the angular field's does — the ONE-space principle of
  ruling R-6.2c-1 (*the carrier's norm is the field's energy*) applied to
  the tail factor;
* the widened ``scalar_flux`` self-derive is honest (the product's cell-group
  factor IS the carrier's widened bulk), so its S4 refusal retired;
* ``SpatialMomentSpace`` is gone; its "append iff > 1" helper is the space
  composition rule :func:`orpheus.numerics.space.factor_tail`.

The witness that ``SpatialMomentSpace`` is unspellable is a numerics claim, so it
lives in ``tests/gates/numerics/test_moment_layout_retired.py``
(``test_the_retired_class_is_unspellable``); this file moved from
``tests/gates/numerics/`` to ``tests/gates/transport/spatial/`` with the layer-DAG
re-filing (2026-10-25), its rows unchanged.

Foundation mark: software invariants of the space layer; no physics claim.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.axis import BasisKind
from orpheus.numerics.basis import SPATIAL_MOMENT_AXIS_LABEL
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace, TensorProductSpace, factor_tail
from orpheus.numerics.coefficient_layout import HarmonicLayout
from orpheus.numerics.spaces.moment_head import head_layout
from orpheus.sn.problem import SNProblem
from orpheus.transport.fields import HarmonicMomentFlux
from orpheus.transport.fields._bases import BulkField
from orpheus.transport.frames import HarmonicFrame
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from tests.gates.sn._test_helpers import placeholder_materials
from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Modal
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


pytestmark = pytest.mark.foundation


def _ld_2d() -> SNProblem:
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4), edges_y=np.linspace(0.0, 1.0, 3),
        mat_map=np.zeros((3, 2), dtype=int), coord=CoordSystem.CARTESIAN,
        face_roles={
            "xmin": MirrorDeck(),
            "xmax": MirrorDeck(),
            "ymin": MirrorDeck(),
            "ymax": MirrorDeck(),
        },
    )
    return SNProblem(mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=2), scheme=LinearDiscontinuous(), boundary_conditions=conditions_for(mesh, _CONDITIONS))


def _ld_1d() -> SNProblem:
    geometry = StructuredGeometry.slab((0.0, 1.0), (0,), left=0, right=0)
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(5)).mesh
    return SNProblem(mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=2), scheme=LinearDiscontinuous(), boundary_conditions=conditions_for(mesh, _CONDITIONS))


def _dd_2d() -> SNProblem:
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4), edges_y=np.linspace(0.0, 1.0, 3),
        mat_map=np.zeros((3, 2), dtype=int), coord=CoordSystem.CARTESIAN,
        face_roles={
            "xmin": MirrorDeck(),
            "xmax": MirrorDeck(),
            "ymin": MirrorDeck(),
            "ymax": MirrorDeck(),
        },
    )
    return SNProblem(mesh, Quadrature.level_symmetric(4), placeholder_materials(ng=2), scheme=DiamondDifference(), boundary_conditions=conditions_for(mesh, _CONDITIONS))


_LD = {"ld_2d": _ld_2d, "ld_1d": _ld_1d}


def _tail_axes(space: FunctionSpace):
    assert space.axes is not None, "a tailed space is axis-built"
    return [ax for ax in space.axes if ax.label == SPATIAL_MOMENT_AXIS_LABEL]


@pytest.mark.parametrize("label", list(_LD), ids=list(_LD))
@pytest.mark.parametrize("L", [0, 1])
def test_the_widened_moment_product_is_axis_built_with_the_schemes_axis_as_its_tail(label: str, L: int) -> None:
    sn = _LD[label]()
    per_axis = sn.scheme.spatial_basis_per_axis
    assert per_axis == 2
    space = sn.moment_space(L, spatial_moments=per_axis)
    assert isinstance(space, TensorProductSpace)
    assert space.axes is not None and space.inner_product_weights is None
    assert space.shape[-1] == per_axis ** sn.ndim
    (tail,) = _tail_axes(space)
    expected = sn.scheme.moment_axis(sn.axes)
    assert tail == expected, "the tail IS the scheme's own moment axis"
    assert isinstance(tail.kind, Modal)
    assert tail.weights is not None
    # The mass value is pinned in
    # test_cell_basis_head.py::test_the_moment_axis_metric_is_the_gram_not_its_inverse.
    # the cell-group factor is the carrier's WIDENED bulk — the same space the scalar family mints
    assert space.factors[1] == BulkField.compose_spatial_moments(sn.bulk_space, sn, per_axis)
    assert space.factors[0] == frame_faces.basis_space(frame_faces.angular_frame(sn.quad, L))
    # width reads
    assert BulkField.spatial_moments_per_axis_of(space) == per_axis
    assert BulkField._spatial_moment_tail_of(space) == (per_axis ** sn.ndim,)
    # one object per (L, width) within the owner
    assert sn.moment_space(L, spatial_moments=per_axis) is space


@pytest.mark.parametrize("label", list(_LD), ids=list(_LD))
def test_the_hub_and_the_frame_agree_at_width_two_as_at_width_one(label: str) -> None:
    sn = _LD[label]()
    L = 1
    widened_angular = BulkField.compose_spatial_moments(sn.angular_bulk_space, sn, 2)
    assert _tail_axes(widened_angular), "the widened angular space carries the scheme's axis"
    frame = HarmonicFrame.for_space(widened_angular, L)
    derived = frame.moment_space_on(widened_angular)
    hub = sn.moment_space(L, spatial_moments=2)
    assert hub == derived and hash(hub) == hash(derived), "two owners, ONE widened space"
    # the tail axis the frame threads IS the angular space's own axis object
    assert derived.axes is not None
    assert _tail_axes(derived)[0] is _tail_axes(widened_angular)[0]
    # and the analysis face's codomain is that space
    assert frame.flux_analysis_on(widened_angular).codomain == hub
    # the un-widened pair still agrees
    assert sn.moment_space(L) == frame.moment_space_on(sn.angular_bulk_space)


@pytest.mark.parametrize("label", list(_LD), ids=list(_LD))
def test_the_moment_norm_carries_the_cell_mass_on_its_tail(label: str) -> None:
    """Parseval across the tail: on a widened iterate the moment space's pairing
    factorises as (head) × (bulk) × (mass-weighted tail) — the tail's measure
    is the scheme's mass, exactly the angular space's. NEGATIVE CONTROL: a
    Euclidean tail (the retired spelling) differs by the mass ratio."""
    sn = _LD[label]()
    space = sn.moment_space(1, spatial_moments=2)
    assert isinstance(space, TensorProductSpace)
    (tail,) = _tail_axes(space)
    assert tail.weights is not None
    rng = np.random.default_rng(3)
    x = rng.standard_normal(space.shape)
    got = space.inner_product(x, x)
    # oracle: the un-widened product's pairing per tail slot, weighted by the mass
    narrow = sn.moment_space(1)
    per_slot = np.array([narrow.inner_product(x[..., k], x[..., k]) for k in range(space.shape[-1])])
    np.testing.assert_allclose(got, float(per_slot @ tail.weights), rtol=1e-12)
    euclidean = float(per_slot.sum())
    assert not np.isclose(got, euclidean), "the mass must move the pairing, or the tail carries none"


def test_width_one_appends_nothing_and_the_tail_rule_is_space_composition() -> None:
    sn = _dd_2d()
    space = sn.moment_space(1)
    assert _tail_axes(space) == []
    assert BulkField.spatial_moments_per_axis_of(space) == 1
    assert sn.moment_space(1, spatial_moments=1) is space
    assert factor_tail(1) == () and factor_tail(4) == (4,)
    # a DD carrier mints no moment axis, so a widened request is refused by the scheme
    with pytest.raises(NotImplementedError, match="no moment axis"):
        sn.moment_space(1, spatial_moments=2)
    # the sphere rule's head is the harmonic family (since S7 read off the
    # head axis's layout; until then ``find_factor`` on the family class),
    # and no other factor of the product is a moment head (the door refuses).
    assert isinstance(space, TensorProductSpace)
    assert isinstance(head_layout(space.factors[0]), HarmonicLayout)
    with pytest.raises(TypeError, match="is not a moment head"):
        head_layout(space.factors[1])


@pytest.mark.parametrize("label", list(_LD), ids=list(_LD))
def test_a_widened_moment_field_self_derives_its_scalar_flux_and_truncates_on_the_hubs_space(label: str) -> None:
    sn = _LD[label]()
    L = 1
    field = HarmonicMomentFlux.zeros_for_problem_and_L(sn, L, spatial_moments=2)
    assert field.space is sn.moment_space(L, spatial_moments=2)
    scalar = field.scalar_flux()
    assert scalar.space == BulkField.compose_spatial_moments(sn.bulk_space, sn, 2)
    assert _tail_axes(scalar.space) and _tail_axes(scalar.space)[0] == sn.scheme.moment_axis(sn.axes)
    truncated = field.truncate(0)
    assert truncated.space == sn.moment_space(0, spatial_moments=2)
    assert truncated.spatial_moments == 2

