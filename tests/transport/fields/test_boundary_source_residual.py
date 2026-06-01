r"""Foundation tests for the B.3 boundary role leaves
:class:`~orpheus.transport.sources.boundary_source.BoundarySource` and
:class:`~orpheus.transport.residuals.boundary_residual.BoundaryResidual`.

These complete the ``{Boundary}×{Source,Residual}`` cells of the field
role grid (siblings of
:class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`), giving
the :class:`~orpheus.transport.fields._bases.BoundaryField` storage base
its 2nd and 3rd concrete instances. All three boundary leaves are empty
leaves — storage / validation / algebra / per-face access / factories
are inherited from ``BoundaryField`` — so most of their machinery is
already pinned by ``test_boundary_flux.py``. This module adds:

* construction of the two NEW leaves (``zeros_for_sn_mesh`` /
  ``from_face_arrays``), and
* the **load-bearing cross-class invariant** unique to the boundary
  family: all three leaves share the SAME ``TraceSpace`` (``mesh.trace``),
  so the space gate would PASS — it is the **class-identity** gate that
  rejects ``BoundarySource ± BoundaryFlux`` etc. This is the boundary
  analogue of the bulk "same units / same space ≠ same meaning"
  discipline.

``foundation`` tests — software invariants, not physics claims.

References
----------

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.3.
* ``orpheus/numerics/field.py`` — Layer-1 class-identity gate.
* ``tests/transport/fields/test_boundary_flux.py`` — the BoundaryFlux
  contract these leaves inherit.
"""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.field import Field
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields._bases import BoundaryField
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.residuals import BoundaryResidual
from orpheus.transport.sources import BoundarySource

from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]

# The two NEW boundary leaves (parametrized over for shared contract).
NEW_BOUNDARY_LEAVES = [BoundarySource, BoundaryResidual]


# ── Fixtures ─────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 5, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere_mesh(nx: int = 5, ng: int = 2) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _cartesian_2d_mesh(nx: int = 3, ny: int = 2, ng: int = 2) -> SNMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ════════════════════════════════════════════════════════════════════
# Construction + inheritance (both new leaves share BoundaryFlux's
# inherited machinery; spot-check it is wired on the new classes).
# ════════════════════════════════════════════════════════════════════


@pytest.mark.parametrize("Leaf", NEW_BOUNDARY_LEAVES)
class TestConstructionInherited:
    def test_inherits_field_and_boundary_field(self, Leaf) -> None:
        m = _slab_mesh()
        bf = Leaf.zeros_for_sn_mesh(m)
        assert isinstance(bf, Field)
        assert isinstance(bf, BoundaryField)
        assert isinstance(bf, Leaf)
        assert bf.mesh is m

    def test_zeros_for_sn_mesh_uses_mesh_trace(self, Leaf) -> None:
        m = _slab_mesh()
        bf = Leaf.zeros_for_sn_mesh(m)
        assert bf.space is m.trace
        np.testing.assert_array_equal(bf.values, 0.0)
        assert set(bf.layout.faces) == {"xmin", "xmax"}

    def test_sphere_layout_only_xmax(self, Leaf) -> None:
        m = _sphere_mesh()
        bf = Leaf.zeros_for_sn_mesh(m)
        assert set(bf.layout.faces) == {"xmax"}

    def test_from_face_arrays_slab(self, Leaf) -> None:
        m = _slab_mesh()
        N = m.quad.N
        bf = Leaf.from_face_arrays(
            m, {"xmin": np.full((N, 2), 1.0), "xmax": np.full((N, 2), 2.0)},
        )
        assert isinstance(bf, Leaf)
        np.testing.assert_array_equal(bf.face_view("xmin"), 1.0)
        np.testing.assert_array_equal(bf.face_view("xmax"), 2.0)

    def test_from_face_arrays_2d(self, Leaf) -> None:
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        N = m.quad.N
        bf = Leaf.from_face_arrays(m, {
            "xmin": np.full((N, 2, 2), 1.0),
            "xmax": np.full((N, 2, 2), 2.0),
            "ymin": np.full((N, 2, 3), 3.0),
            "ymax": np.full((N, 2, 3), 4.0),
        })
        np.testing.assert_array_equal(bf.face_view("ymin"), 3.0)

    def test_post_init_requires_trace_space(self, Leaf) -> None:
        m = _slab_mesh()
        from orpheus.numerics.space import FunctionSpace
        plain = FunctionSpace(name="sn_boundary_flat", shape=(m.trace.shape[0],))
        with pytest.raises(TypeError, match="TraceSpace"):
            Leaf(values=np.zeros(m.trace.shape[0]), space=plain, mesh=m)


@pytest.mark.parametrize("Leaf", NEW_BOUNDARY_LEAVES)
class TestAlgebraClosedWithinClass:
    def test_add_sub_closed(self, Leaf) -> None:
        m = _slab_mesh()
        a = Leaf.zeros_for_sn_mesh(m)
        b = Leaf.zeros_for_sn_mesh(m)
        a.values[:] = 1.0
        b.values[:] = 2.0
        s = a + b
        d = b - a
        assert isinstance(s, Leaf) and isinstance(d, Leaf)
        np.testing.assert_array_equal(s.values, 3.0)
        np.testing.assert_array_equal(d.values, 1.0)

    def test_scalar_mul_div_neg(self, Leaf) -> None:
        m = _sphere_mesh()
        a = Leaf.zeros_for_sn_mesh(m)
        a.values[:] = 4.0
        np.testing.assert_allclose((2.5 * a).values, 10.0)
        np.testing.assert_allclose((a / 4.0).values, 1.0)
        np.testing.assert_array_equal((-a).values, -4.0)
        assert isinstance(-a, Leaf)

    def test_cross_mesh_add_rejected(self, Leaf) -> None:
        a = Leaf.zeros_for_sn_mesh(_slab_mesh())
        b = Leaf.zeros_for_sn_mesh(_slab_mesh())  # distinct instance
        with pytest.raises(ValueError, match="mesh-bound"):
            _ = a + b

    def test_frozen(self, Leaf) -> None:
        m = _slab_mesh()
        bf = Leaf.zeros_for_sn_mesh(m)
        with pytest.raises(FrozenInstanceError):
            bf.mesh = m  # type: ignore[misc]


# ════════════════════════════════════════════════════════════════════
# The load-bearing boundary invariant: cross-class arithmetic RAISES
# even though all three boundary leaves share the SAME TraceSpace.
#
# This is structurally DIFFERENT from the bulk families, where
# cross-class operands also have different spaces/shapes (so either gate
# would catch them). Here the space gate would PASS — only the
# class-identity gate discriminates BoundaryFlux / BoundarySource /
# BoundaryResidual.
# ════════════════════════════════════════════════════════════════════


class TestCrossClassRejectionSharedSpace:
    def test_all_three_share_the_same_trace_space(self) -> None:
        """Pre-condition for the invariant: the three leaves are built on
        the IDENTICAL ``mesh.trace`` object — so ``space == space`` and
        the space gate alone would NOT reject cross-class arithmetic."""
        m = _slab_mesh()
        flux = BoundaryFlux.zeros_for_sn_mesh(m)
        src = BoundarySource.zeros_for_sn_mesh(m)
        res = BoundaryResidual.zeros_for_sn_mesh(m)
        assert flux.space is src.space is res.space is m.trace
        # The space gate would pass (equal spaces) — proving it is the
        # CLASS gate that must do the rejection below.
        assert flux.space == src.space == res.space

    @pytest.mark.parametrize(
        ("A", "B"),
        [
            (BoundarySource, BoundaryFlux),
            (BoundaryResidual, BoundaryFlux),
            (BoundarySource, BoundaryResidual),
        ],
    )
    def test_cross_class_add_sub_raises(self, A, B) -> None:
        m = _slab_mesh()
        a = A.zeros_for_sn_mesh(m)
        b = B.zeros_for_sn_mesh(m)
        # Same shape, same space, distinct class → class gate rejects.
        assert a.values.shape == b.values.shape
        assert a.space == b.space
        with pytest.raises(TypeError, match="same-class"):
            _ = a + b  # type: ignore[operator]
        with pytest.raises(TypeError, match="same-class"):
            _ = b - a  # type: ignore[operator]

    def test_cross_class_inner_product_raises(self) -> None:
        m = _slab_mesh()
        src = BoundarySource.zeros_for_sn_mesh(m)
        res = BoundaryResidual.zeros_for_sn_mesh(m)
        with pytest.raises(TypeError, match="same-class"):
            src.inner_product(res)  # type: ignore[arg-type]
