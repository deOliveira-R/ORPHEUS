r"""One moment space, two owners — the hub and the frame agree STRUCTURALLY
(CS4c step 6 item 6.2c-ii; the memo's gates P5, P7 and P10).

Since item 6.2b the carrier (:class:`~orpheus.sn.problem.SNProblem`)
owns the moment space and every moment field / admission guard on it holds
ONE object per ``(L, width)``; since item 6.2c-ii the head it composes is the
frame's Parseval-dressed ``basis_space`` (ruling R-6.2c-1), and the frame's
own derivation, :meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.moment_space_on`,
builds the same product from the angular space's own axes. The two owners
cannot share an object (`[M]` ruling O-5: the frame has no carrier to reach),
so the reachable target — and the one asserted here — is structural equality
ACROSS owners and ``is``-identity WITHIN each.

* **P5** — hub ``==`` frame, ``is`` within each, on every ledger carrier;
* **P7** — the three product mechanisms are ONE space by identity — the
  hub's ``head * bulk``, ``FunctionSpace.of_axes(head_axis, *bulk.axes)``,
  and the frame's derivation — while their NAMES differ (`[M]` the ``*``
  product keeps the head's family tag, the ``of_axes`` mint is a pure
  digest), so the gate asserts identity and NOT the name;
* **P10** — the structural / memory leg: on every (geometry × L) row the
  moment space is axis-built, ``inner_product_weights is None``, the metric
  is derived (no object on a DIAGONAL-Gram frame), and no reachable array
  has the product's size.

Foundation mark: software invariants (identity, structure); no physics claim.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.frame import GramStructure
from orpheus.numerics.metric import FactoredMetric
from orpheus.numerics.space import FunctionSpace, TensorProductSpace
from orpheus.transport.frames import HarmonicFrame
from tests.gates.sn.architecture.test_monomorphic_leaves import (
    _cart2d,
    _cylinder,
    _slab,
    _sphere,
)
from orpheus.numerics import frame_faces
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.basis import SPATIAL_MOMENT_AXIS_LABEL
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.fields._bases import BulkField
from orpheus.transport.spatial import LinearDiscontinuous
from tests.gates.transport._materials import placeholder_materials
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

pytestmark = pytest.mark.foundation

_GEOMETRIES = {"slab": _slab, "sphere": _sphere, "cylinder": _cylinder, "cart2d": _cart2d}
_ORDERS = [0, 1, 2]


@pytest.mark.parametrize("geometry", list(_GEOMETRIES), ids=list(_GEOMETRIES))
@pytest.mark.parametrize("L", _ORDERS)
def test_P5_the_hub_and_the_frame_derive_one_space_by_equality_and_each_owner_by_identity(
    geometry: str, L: int,
) -> None:
    problem = _GEOMETRIES[geometry]()
    hub = problem.moment_space(L)
    frame = HarmonicFrame.for_space(problem.angular_bulk_space, L)
    derived = frame.moment_space_on(problem.angular_bulk_space)
    assert hub == derived and derived == hub, "two owners, ONE space (structural equality)"
    assert hash(hub) == hash(derived)
    # is-identity WITHIN each owner
    assert problem.moment_space(L) is hub
    assert frame.flux_analysis_on(problem.angular_bulk_space).codomain == hub
    # the head both owners hold is the frame's Parseval-dressed one
    assert isinstance(hub, TensorProductSpace)
    assert hub.factors[0] == frame_faces.basis_space(frame)
    assert hub.factors[0] != FunctionSpace.of_basis(frame.basis), "the continuum head is another space"


@pytest.mark.parametrize("geometry", list(_GEOMETRIES), ids=list(_GEOMETRIES))
@pytest.mark.parametrize("L", _ORDERS)
def test_P7_the_three_product_mechanisms_are_one_space_by_identity_not_by_name(
    geometry: str, L: int,
) -> None:
    problem = _GEOMETRIES[geometry]()
    hub = problem.moment_space(L)
    frame = frame_faces.angular_frame(problem.quad, L)
    head = frame_faces.basis_space(frame)
    bulk = problem.bulk_space
    assert head.axes is not None and bulk.axes is not None
    star = head * bulk
    of_axes = FunctionSpace.of_axes(*head.axes, *bulk.axes)
    derived = HarmonicFrame.from_galerkin(frame).moment_space_on(problem.angular_bulk_space)
    assert hub == star == of_axes == derived
    assert len({hub, star, of_axes, derived}) == 1
    # the names DIFFER — identity is structural, never nominal
    assert star.name == hub.name
    assert of_axes.name != hub.name, "the of_axes mint is a pure digest; the * product keeps the head's tag"
    # the DENSE-Gram overlay (if any) rides both mechanisms alike
    if frame.discrete_gram_structure is GramStructure.DENSE:
        assert isinstance(star.metric, FactoredMetric)
    else:
        assert star.metric is None
    x = np.random.default_rng(L).standard_normal(hub.shape)
    np.testing.assert_array_equal(star.apply_metric(x), of_axes.apply_metric(x))
    np.testing.assert_array_equal(star.apply_metric(x), derived.apply_metric(x))


@pytest.mark.parametrize("geometry", list(_GEOMETRIES), ids=list(_GEOMETRIES))
@pytest.mark.parametrize("L", _ORDERS)
def test_P10_the_moment_space_is_axis_built_and_never_densified(geometry: str, L: int) -> None:
    problem = _GEOMETRIES[geometry]()
    space = problem.moment_space(L)
    frame = frame_faces.angular_frame(problem.quad, L)
    assert space.axes is not None, "the moment space is axis-built"
    assert space.inner_product_weights is None
    if frame.discrete_gram_structure is GramStructure.DENSE:
        assert isinstance(space.metric, FactoredMetric)
        forms = [f for _, f in space.metric.entries]
        assert forms[0] is not None and all(f is None for f in forms[1:])
    else:
        assert space.metric is None
    n = int(np.prod(space.shape))
    reachable = [ax.weights for ax in space.axes if ax.weights is not None]
    assert all(w.size < n for w in reachable), "no reachable array is state-sized"
    assert sum(w.nbytes for w in reachable) < 8 * n


# ── Moved from tests/gates/transport/spatial/test_spatial_moment_tail_is_the_schemes_axis.py by layer DAG Phase 4
# (2026-11-01): the hub and the frame agree at width two as at width one (the hub's keyed moment space against the
# frame's mint of the widened angular space), with its two LD problems.


#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


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


_LD = {"ld_2d": _ld_2d, "ld_1d": _ld_1d}


def _tail_axes(space: FunctionSpace):
    assert space.axes is not None, "a tailed space is axis-built"
    return [ax for ax in space.axes if ax.label == SPATIAL_MOMENT_AXIS_LABEL]


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
