r"""``scheme.trial_space(space, axes)``: a space in the closure's within-cell basis (layer DAG Phase 4).

The verb is the one spelling of "append the scheme's moment axis". Its law, posed on axis-built spaces built here
from numerics axes and a mesh's axis tuple:

* a slopeless closure (Diamond Difference) returns ``space`` ITSELF, the same object (`is`), on every chart and on
  an axes-less space too: its trial basis is the cell average;
* a multi-moment closure (Linear Discontinuous) returns ``space.with_axis(moment_axis(axes))``: the leading
  axes are the input's own axis objects, and the appended axis carries the scheme's mass, the hand reference
  ``prod_a 1/(2 o_a + 1)`` over the Legendre orders ``o_a in {0, 1}``: ``[1, 1/3]`` on a slab and
  ``[1, 1/3, 1/3, 1/9]`` on a 2-D Cartesian mesh. The LD row's result is a fresh value: ``is`` is NOT asserted;
  equality to a second call is (a value type re-mints);
* Linear Discontinuous refuses an axes-less space (``TypeError`` from ``FunctionSpace.with_axis``, fragment "is not axis-built") and a radial
  chart (``NotImplementedError`` from ``moment_axis``, fragment "has no moment mass on a radial_"; that refusal is
  the scheme's, pinned at its source by ``test_cell_basis_head.py``; this row pins that the verb does not
  swallow it).

The population is DD and LD only: ``Step`` is named in the verb's docstring but has no class in the tree
(``[M]`` 2026-11-01, the census's ``__subclasses__`` walk). Named first reds: a verb that always widens (the DD
rows), a verb that copies on a slopeless closure (the DD ``is`` rows), and a Euclidean moment axis (the LD weight
rows).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.geometry.deck import MirrorDeck
from orpheus.mesh import CellsByCount, Mesh2D, Mesher, axes_from_legacy_mesh
from orpheus.numerics.axis import Axis, BasisKind, Declared
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous

pytestmark = pytest.mark.foundation


def _slab_axes():
    g = StructuredGeometry.slab((0.0, 1.0), (0,), left=MirrorDeck(), right=MirrorDeck())
    return axes_from_legacy_mesh(Mesher(g).partition(CellsByCount.uniform_width(4)).mesh)


def _cart2d_axes():
    return axes_from_legacy_mesh(Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4), edges_y=np.linspace(0.0, 1.0, 3), mat_map=np.zeros((3, 2), dtype=int),
        face_roles={face: MirrorDeck() for face in ("xmin", "xmax", "ymin", "ymax")},
    ))


def _sphere_axes():
    g = StructuredGeometry.sphere((0.0, 1.0), (0,), outer=0)
    return axes_from_legacy_mesh(Mesher(g).partition(CellsByCount.uniform_width(4)).mesh)


_AXES = {"slab": _slab_axes, "cart2d": _cart2d_axes, "sphere": _sphere_axes}


def _space(axes) -> FunctionSpace:
    """An axis-built angular bulk: ordinates x two groups x the spatial cells, with a non-uniform spatial weight."""
    spatial = tuple(a.n for a in axes)
    w = np.linspace(0.5, 2.0, int(np.prod(spatial))).reshape(spatial)
    return FunctionSpace.of_axes(
        Axis.of_quadrature(Quadrature.gauss_legendre(4)),
        Axis("energy", (2,), origin=Declared(BasisKind.NODAL)),
        Axis("spatial", spatial, weights=w, origin=Declared(BasisKind.NODAL)),
    )


#: The hand mass of the LD moment axis: prod over axes of 1/(2 o + 1), o in {0, 1}, slot order of the tensor basis.
_LD_MASS = {"slab": np.array([1.0, 1.0 / 3.0]), "cart2d": np.array([1.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 9.0])}


@pytest.mark.parametrize("chart", list(_AXES))
def test_a_slopeless_closure_returns_the_space_itself(chart: str) -> None:
    axes = _AXES[chart]()
    space = _space(axes)
    if DiamondDifference().trial_space(space, axes) is not space:
        pytest.fail(f"{chart}: DD's trial space is not its input object")


def test_a_slopeless_closure_returns_an_axes_less_space_itself() -> None:
    bare = FunctionSpace("p4_bare", (4, 2, 4))
    if DiamondDifference().trial_space(bare, _slab_axes()) is not bare:
        pytest.fail("DD's trial space of an axes-less space is not its input object")


@pytest.mark.parametrize("chart", list(_LD_MASS))
def test_a_multi_moment_closure_appends_the_moment_axis_with_its_mass(chart: str) -> None:
    axes = _AXES[chart]()
    space = _space(axes)
    trial = LinearDiscontinuous().trial_space(space, axes)
    trial_axes, space_axes = trial.axes, space.axes
    if trial_axes is None or space_axes is None or len(trial_axes) != len(space_axes) + 1:
        pytest.fail(f"{chart}: LD's trial space is not the input's axes plus one")
    for k, (a, b) in enumerate(zip(trial_axes, space_axes)):
        if a is not b:
            pytest.fail(f"{chart}: leading axis {k} is not the input's own axis object")
    tail = trial_axes[-1]
    if tail.label != "spatial_moment":
        pytest.fail(f"{chart}: the appended axis is labelled {tail.label!r}")
    np.testing.assert_allclose(np.asarray(tail.weights, float).ravel(), _LD_MASS[chart], rtol=0, atol=4 * np.finfo(float).eps)
    if trial.shape != (*space.shape, _LD_MASS[chart].size):
        pytest.fail(f"{chart}: trial shape {trial.shape}")
    again = LinearDiscontinuous().trial_space(space, axes)
    if again != trial or hash(again) != hash(trial):
        pytest.fail(f"{chart}: two calls give unequal trial spaces")


def test_a_multi_moment_closure_refuses_an_axes_less_space() -> None:
    with pytest.raises(TypeError, match="is not axis-built"):
        LinearDiscontinuous().trial_space(FunctionSpace("p4_bare", (4, 2, 4)), _slab_axes())


@pytest.mark.rests_on("tests/gates/transport/spatial/test_cell_basis_head.py::test_the_ld_moment_axis_is_the_cell_basis_head")
def test_a_multi_moment_closure_on_a_radial_chart_refuses_through_the_verb() -> None:
    axes = _sphere_axes()
    with pytest.raises(NotImplementedError, match="has no moment mass on a radial_"):
        LinearDiscontinuous().trial_space(_space(axes), axes)


@pytest.mark.parametrize(("space_chart", "axes_chart"), [("slab", "cart2d"), ("cart2d", "slab")], ids=["slab-space-cart2d-axes", "cart2d-space-slab-axes"])
def test_a_multi_moment_closure_refuses_axes_that_did_not_build_the_space(space_chart: str, axes_chart: str) -> None:
    """The moment axis reads the mesh axes' count, so axes other than those the space's spatial axis was built over
    would widen the space with another dimension's cell basis: LD refuses the pair by name. Both directions (a slab
    bulk with rectangle axes, a rectangle bulk with slab axes); a slopeless closure has nothing to widen and returns
    the space (the rows above)."""
    with pytest.raises(ValueError, match="trial_space: the space's spatial axis has shape"):
        LinearDiscontinuous().trial_space(_space(_AXES[space_chart]()), _AXES[axes_chart]())
