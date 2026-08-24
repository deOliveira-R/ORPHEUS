r"""CS4b S4 — fields are space elements: the mesh binding is GONE (G4.1/G4.2).

The deletion step's own witnesses (verification plan §10 S4; done-when #1
and #2):

* **G4.1 (the headline)** — a leaf constructs from ``(values, space)``
  alone. ``[M]`` before S4 this exact call raised ``TypeError: missing 1
  required keyword-only argument: 'mesh'``.
* **G4.2** — the ``mesh`` dataclass field is absent from every shipped
  role leaf (the two declaration roots and their six covariant
  narrowings died together), and the composite exposes no ``mesh``
  surface — the carrier's knowledge enters through its cached space
  mints (read at the call site, S5) and the ``space_on`` admission
  seams only.

* **G5.1 (the S5 retirement gate)** — the mesh-keyed sugar tier
  (``from_mesh`` / ``zeros_on`` / ``from_ndarray``) is GONE from every
  concrete field leaf, and the composite allocators are space-keyed.
  The named survivors are the NON-field tiers: ``MaterialXSField.
  from_mesh`` (mesh→container assembly), ``MultiplicationOperator.
  from_mesh`` (operator tier), the BC factories, and the S6-pending
  moment-family keyed factories.

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


class TestG51SugarTierRetired:
    """G5.1 — the S5 done-when as a permanent gate.

    The step's grep predicate ("no leaf spells the mesh-keyed sugar"),
    pinned structurally so a nostalgic re-addition reds here: the sugar
    names are ABSENT from every concrete leaf, the composite allocators
    are space-keyed, and the documented survivors still exist (so this
    gate cannot silently widen into banning the non-field tiers).
    """

    _LEAVES = ()  # populated below — import-time failures stay readable

    def _all_leaves(self):
        from orpheus.transport.fields.angular_flux import AngularFlux
        from orpheus.transport.fields.scalar_flux import ScalarFlux
        from orpheus.transport.fields.cross_section_field import (
            CrossSectionField,
        )
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )
        from orpheus.transport.fields.scalar_boundary_flux import (
            ScalarBoundaryFlux,
        )
        from orpheus.transport.fields.harmonic_moment_flux import (
            HarmonicMomentFlux,
        )
        from orpheus.transport.fields.radial_characteristic_interior_flux import (
            RadialCharacteristicInteriorFlux,
        )
        from orpheus.transport.fields.radial_characteristic_boundary_flux import (
            RadialCharacteristicBoundaryFlux,
        )
        from orpheus.transport.source_sinks import (
            AngularBoundarySourceSink,
            AngularSourceSink,
            ScalarBoundarySourceSink,
            ScalarSourceSink,
        )
        from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
            RadialCharacteristicInteriorSourceSink,
        )
        from orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink import (
            RadialCharacteristicBoundarySourceSink,
        )
        from orpheus.transport.residuals import (
            AngularBoundaryResidual,
            AngularResidual,
            ScalarResidual,
        )
        from orpheus.transport.residuals.radial_characteristic_interior_residual import (
            RadialCharacteristicInteriorResidual,
        )
        from orpheus.transport.residuals.radial_characteristic_boundary_residual import (
            RadialCharacteristicBoundaryResidual,
        )

        return (
            AngularFlux, ScalarFlux, CrossSectionField,
            AngularSourceSink, ScalarSourceSink,
            AngularResidual, ScalarResidual,
            AngularBoundaryFlux, AngularBoundarySourceSink,
            AngularBoundaryResidual,
            ScalarBoundaryFlux, ScalarBoundarySourceSink,
            HarmonicMomentFlux,
            RadialCharacteristicInteriorFlux,
            RadialCharacteristicInteriorSourceSink,
            RadialCharacteristicInteriorResidual,
            RadialCharacteristicBoundaryFlux,
            RadialCharacteristicBoundarySourceSink,
            RadialCharacteristicBoundaryResidual,
        )

    def test_no_concrete_leaf_exposes_the_sugar_tier(self) -> None:
        offenders = [
            f"{leaf.__name__}.{name}"
            for leaf in self._all_leaves()
            for name in ("from_mesh", "zeros_on", "from_ndarray")
            if hasattr(leaf, name)
        ]
        if offenders:
            pytest.fail(
                "the S5-retired sugar tier resurfaced on: "
                + ", ".join(offenders)
            )

    def test_composite_allocators_are_space_keyed(self) -> None:
        import inspect

        from orpheus.transport.full_field import FullField
        from orpheus.transport.radial_characteristic_field import (
            RadialCharacteristicField,
        )
        from orpheus.transport.timed_full_field import TimedFullField

        for func in (FullField.zeros, TimedFullField.zeros):
            params = inspect.signature(func).parameters
            if "space" not in params or "mesh" in params:
                pytest.fail(f"{func.__qualname__} is not space-keyed")
        for func in (
            RadialCharacteristicField.flux_zeros,
            RadialCharacteristicField.source_zeros,
        ):
            params = inspect.signature(func).parameters
            if "space" not in params:
                pytest.fail(f"{func.__qualname__} is not space-keyed")

    def test_the_documented_survivors_still_exist(self) -> None:
        """The non-field tiers this gate deliberately does NOT ban."""
        from orpheus.transport.mesh.material_xs_field import MaterialXSField
        from orpheus.transport.operators.multiplication_operator import (
            MultiplicationOperator,
        )
        from orpheus.transport.fields.harmonic_moment_flux import (
            HarmonicMomentFlux,
        )

        if not callable(getattr(MaterialXSField, "from_mesh", None)):
            pytest.fail("MaterialXSField.from_mesh (assembly tier) missing")
        if not callable(getattr(MultiplicationOperator, "from_mesh", None)):
            pytest.fail("MultiplicationOperator.from_mesh (operator tier) missing")
        if not callable(getattr(HarmonicMomentFlux, "from_mesh_and_L", None)):
            pytest.fail("the moment family's keyed factory missing (S6 scope)")
