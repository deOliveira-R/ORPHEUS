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

Layer DAG Phase 4 (2026-11-01): only the Parseval-across-the-tail row stays here, posed on material meshes with the
frame's moment space; the rows on the SN problem's keyed moment space moved to
tests/gates/sn/mesh/test_spatial_moment_tail_on_the_problem.py and test_hub_and_frame_agree_on_the_moment_space.py.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.basis import SPATIAL_MOMENT_AXIS_LABEL
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace, TensorProductSpace
from orpheus.transport.mesh.material_mesh import MaterialMesh
from tests.gates.transport._phase_space_cases import Posed, moment_space
from orpheus.transport.spatial import LinearDiscontinuous
from tests.gates.transport._materials import placeholder_materials
from orpheus.geometry.deck import MirrorDeck


pytestmark = pytest.mark.foundation


def _ld_2d() -> Posed:
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
    return Posed(MaterialMesh(mesh, placeholder_materials(ng=2)), Quadrature.level_symmetric(4), LinearDiscontinuous())


def _ld_1d() -> Posed:
    geometry = StructuredGeometry.slab((0.0, 1.0), (0,), left=0, right=0)
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(5)).mesh
    return Posed(MaterialMesh(mesh, placeholder_materials(ng=2)), Quadrature.gauss_legendre(4), LinearDiscontinuous())


_LD = {"ld_2d": _ld_2d, "ld_1d": _ld_1d}


def _tail_axes(space: FunctionSpace):
    assert space.axes is not None, "a tailed space is axis-built"
    return [ax for ax in space.axes if ax.label == SPATIAL_MOMENT_AXIS_LABEL]


@pytest.mark.parametrize("label", list(_LD), ids=list(_LD))
def test_the_moment_norm_carries_the_cell_mass_on_its_tail(label: str) -> None:
    """Parseval across the tail: on a widened iterate the moment space's pairing
    factorises as (head) × (bulk) × (mass-weighted tail) — the tail's measure
    is the scheme's mass, exactly the angular space's. NEGATIVE CONTROL: a
    Euclidean tail (the retired spelling) differs by the mass ratio."""
    sn = _LD[label]()
    space = moment_space(sn.mesh, sn.quadrature, 1, sn.scheme)
    assert isinstance(space, TensorProductSpace)
    (tail,) = _tail_axes(space)
    assert tail.weights is not None
    rng = np.random.default_rng(3)
    x = rng.standard_normal(space.shape)
    got = space.inner_product(x, x)
    # oracle: the un-widened product's pairing per tail slot, weighted by the mass
    narrow = moment_space(sn.mesh, sn.quadrature, 1)
    per_slot = np.array([narrow.inner_product(x[..., k], x[..., k]) for k in range(space.shape[-1])])
    np.testing.assert_allclose(got, float(per_slot @ tail.weights), rtol=1e-12)
    euclidean = float(per_slot.sum())
    assert not np.isclose(got, euclidean), "the mass must move the pairing, or the tail carries none"


