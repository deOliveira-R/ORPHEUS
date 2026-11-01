r"""Phase spaces of the shapes production builds, spelled in numerics vocabulary (layer DAG Phase 4).

The numerics law files pose their laws on realistic spaces: an ordinate axis of a shipped quadrature, a synthetic
energy axis and a spatial axis carrying a mesh's cell measure, the linear-discontinuous cell basis as a moment axis,
a frame-dressed harmonic head. Until Phase 4 they read those spaces off ``SNProblem`` carriers, which a numerics gate
may not import; here each is built from numerics (and geometry) objects only. The shapes match production's:
``tests/gates/transport/test_phase_spaces.py::test_the_numerics_twins_have_the_production_shapes`` pins that these
twins have the axis labels, shapes, origin kinds and ordinate axis of ``orpheus.transport.phase_spaces``'s products on
the four monomorphic meshes (spatial weights to 8 eps), and ``tests/gates/sn/mesh/test_hub_and_frame_agree_on_the_moment_space.py`` that the moment space is the
frame head times the cell group. Nothing here asserts equality with a production carrier (a hand volume formula
and the mesh's own arithmetic differ by about 1e-15 relative on a slab, so such a bridge would pin float arithmetic).
"""
from __future__ import annotations

import math
from typing import NamedTuple

import numpy as np

from orpheus.geometry.manifold import RealSpace
from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
from orpheus.numerics.basis.interval_legendre_basis import spatial_moment_basis
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace

#: The literal edges of the monomorphic population (``tests/gates/transport/_phase_space_cases.py``'s
#: ``NONUNIFORM_EDGES`` and two-region rectangle; a numerics file may not import that transport helper).
MONOMORPHIC_EDGES = np.array([0.0, 0.12, 0.35, 0.80, 1.30, 2.00])
MONOMORPHIC_RECTANGLE = (np.array([0.0, 0.1, 0.35, 0.9, 1.6]), np.array([0.0, 0.2, 0.5, 1.4]))


def cell_volumes(chart: str, edges: np.ndarray) -> np.ndarray:
    """Cell volumes of a 1-D chart from its edges: slab ``h``, cylinder ``pi d(r^2)``, sphere ``4/3 pi d(r^3)``."""
    e = np.asarray(edges, dtype=float)
    match chart:
        case "slab":
            return np.diff(e)
        case "cylinder":
            return math.pi * np.diff(e ** 2)
        case "sphere":
            return 4.0 * math.pi / 3.0 * np.diff(e ** 3)
    raise ValueError(chart)


def spatial_axis_1d(chart: str, edges: np.ndarray) -> Axis:
    """The rank-1 spatial axis a 1-D mesh mints: the cell-volume measure on the cell centres, on the real line."""
    e = np.asarray(edges, dtype=float)
    measure = DiscreteMeasure(nodes=0.5 * (e[1:] + e[:-1]), weights=cell_volumes(chart, e), support=RealSpace(1))
    return Axis.of_measure(measure, "spatial")


def spatial_axis_2d(edges_x: np.ndarray, edges_y: np.ndarray) -> Axis:
    """The rank-2 spatial axis a rectangle mints: the cell areas ``dx dy``, declared nodal (rank-d axes carry no
    measure generator)."""
    dx, dy = np.diff(np.asarray(edges_x, float)), np.diff(np.asarray(edges_y, float))
    return Axis("spatial", (dx.size, dy.size), weights=np.outer(dx, dy), origin=Declared(BasisKind.NODAL))


class Twin(NamedTuple):
    """One phase space's factors: the quadrature, the scalar bulk ``energy x spatial`` and the angular bulk."""

    quadrature: Quadrature
    bulk: FunctionSpace
    angular_bulk: FunctionSpace

    def ld_trial(self) -> FunctionSpace:
        """The angular bulk widened by the linear-discontinuous cell basis (2 polynomials per spatial dimension)."""
        ndim = len(self.bulk.axis("spatial").shape)
        return self.angular_bulk.with_axis(Axis.of_basis(spatial_moment_basis(2, ndim)))

    def moment(self, L: int) -> FunctionSpace:
        """The harmonic moment space at order ``L``: the frame-dressed head times the scalar bulk."""
        return frame_faces.basis_space(frame_faces.angular_frame(self.quadrature, L)) * self.bulk


def twin(quadrature: Quadrature, spatial: Axis, ng: int = 2) -> Twin:
    return Twin(
        quadrature,
        FunctionSpace.of_axes(EnergyAxis.synthetic(ng), spatial),
        FunctionSpace.of_axes(Axis.of_quadrature(quadrature), EnergyAxis.synthetic(ng), spatial),
    )
