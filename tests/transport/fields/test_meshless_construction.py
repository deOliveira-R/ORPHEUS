r"""CS4b S4 — fields are space elements: the mesh binding is GONE (G4.1/G4.2).

The deletion step's own witnesses (verification plan §10 S4; done-when #1
and #2):

* **G4.1 (the headline)** — a leaf constructs from ``(values, space)``
  alone. ``[M]`` before S4 this exact call raised ``TypeError: missing 1
  required keyword-only argument: 'mesh'``.
* **G4.2** — the ``mesh`` dataclass field is absent from every shipped
  role leaf (the two declaration roots and their six covariant
  narrowings died together), and the composite exposes no ``mesh``
  surface — the carrier's knowledge enters through factory ARGUMENTS
  and the ``space_on`` admission seams only.

The Q1 flip's witness (the residual rides the operands' shared mint)
lives with the residual tests (``test_typed_residuals`` /
``test_typed_residual_evaluation``); the derived composite space's
witnesses live with the composite tests (the twin-carrier legs).
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.cross_section_field import CrossSectionField
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.residuals import AngularResidual
from orpheus.transport.source_sinks import AngularSourceSink, ScalarSourceSink
from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]


def _slab(nx: int = 4, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=ng))


class TestG41HeadlineConstruction:
    def test_cross_section_field_constructs_from_values_and_space(self) -> None:
        """The step's headline gate: ``CrossSectionField(values=…, space=…)``
        CONSTRUCTS — [M] pre-S4 it raised the missing-``mesh`` TypeError."""
        m = _slab()
        sigma = CrossSectionField(
            values=np.ones(m.bulk_space.shape), space=m.bulk_space,
        )
        assert sigma.space is m.bulk_space
        assert sigma.ng == m.ng

    def test_every_bulk_family_constructs_meshless(self) -> None:
        m = _slab()
        for cls, space in (
            (AngularFlux, m.angular_bulk_space),
            (AngularSourceSink, m.angular_bulk_space),
            (ScalarFlux, m.bulk_space),
            (ScalarSourceSink, m.bulk_space),
            (AngularResidual, m.angular_bulk_space),
        ):
            f = cls(values=np.zeros(space.shape), space=space)
            assert f.space is space, cls.__name__

    def test_face_family_constructs_meshless(self) -> None:
        m = _slab()
        bf = AngularBoundaryFlux(
            values=np.zeros(m.angular_trace.shape), space=m.angular_trace,
        )
        assert bf.space is m.angular_trace
        assert bf.layout is m.angular_trace.layout  # read off the space


class TestG42TheBindingIsGone:
    def test_no_leaf_declares_a_mesh_field(self) -> None:
        """The two declaration roots and their narrowings died together —
        no shipped role leaf carries ``mesh`` in its dataclass fields."""
        for cls in (
            AngularFlux, AngularSourceSink, AngularResidual,
            ScalarFlux, ScalarSourceSink, CrossSectionField,
            AngularBoundaryFlux, HarmonicMomentFlux,
        ):
            names = {f.name for f in dataclasses.fields(cls)}
            assert "mesh" not in names, cls.__name__

    def test_instances_expose_no_mesh_attribute(self) -> None:
        m = _slab()
        psi = AngularFlux.zeros(m.angular_bulk_space)
        assert not hasattr(psi, "mesh")
        comp = FullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, space=m.full_field_space,
        )
        assert not hasattr(comp, "mesh")  # the composite property retired too
