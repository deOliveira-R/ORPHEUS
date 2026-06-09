r"""Foundation tests for the B.3 typed residual role leaves.

Pins the structural contract of
:class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
and
:class:`~orpheus.transport.residuals.scalar_residual.ScalarResidual` —
shape validation, same-class dunder algebra (closed), the mesh-binding
guard, and the **load-bearing cross-class rejection**: a residual and a
source that share the same shape AND the same units still cannot be
added, because Field's Layer-1 class-identity gate makes "same units"
grant permission, not meaning (the named-composition discipline, each
residual leaf's ``from_balance`` factory, B.5 / Issue #201).

These are ``foundation`` tests: they verify software invariants
(constructor shape check, frozen-ness, dunder algebra, the class gate)
rather than physics claims — there is no theory-page ``:label:`` to
verify against.

References
----------

* ``.claude/plans/field_role_typing_view_g.md`` — Phase B step B.3.
* ``orpheus/numerics/field.py`` — Layer-1 class-identity gate
  (``_check_partner``).
"""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.field import Field
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.transport.residuals import (
    AngularResidual,
    BoundaryResidual,
    ScalarResidual,
)
from orpheus.transport.source_sinks import AngularSourceSink, ScalarSourceSink

from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# ── Fixtures ─────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Build a small slab :class:`SNMesh` for unit testing."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _2d_mesh(nx: int = 3, ny: int = 3, ng: int = 1) -> SNMesh:
    """Build a small 2-D Cartesian :class:`SNMesh`."""
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _ang_shape(m: SNMesh) -> tuple[int, ...]:
    return (m.quad.N, m.ng, *m.spatial_shape)


def _sca_shape(m: SNMesh) -> tuple[int, ...]:
    return (m.ng, *m.spatial_shape)


# ════════════════════════════════════════════════════════════════════
# AngularResidual — construction + within-class algebra
# ════════════════════════════════════════════════════════════════════


class TestAngularResidual:
    def test_inherits_field_and_angular_family(self) -> None:
        m = _slab_mesh()
        r = AngularResidual.from_mesh(np.zeros(_ang_shape(m)), m)
        assert isinstance(r, Field)
        assert r.space.name == "angular_residual"
        assert r.mesh is m

    def test_from_mesh_shape_and_metadata(self) -> None:
        m = _slab_mesh()
        r = AngularResidual.from_mesh(np.ones(_ang_shape(m)), m)
        assert r.values.shape == _ang_shape(m)
        assert (r.N, r.ng, r.nx, r.ny) == (m.quad.N, m.ng, m.nx, m.ny)

    def test_from_ndarray_alias(self) -> None:
        m = _slab_mesh()
        r = AngularResidual.from_ndarray(np.ones(_ang_shape(m)), m)
        assert isinstance(r, AngularResidual)

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="AngularResidual.*does not match"):
            AngularResidual.from_mesh(
                np.zeros((m.quad.N + 1, m.ng, *m.spatial_shape)), m,
            )

    def test_within_class_add_sub_closed(self) -> None:
        m = _slab_mesh()
        a = AngularResidual.from_mesh(np.ones(_ang_shape(m)), m)
        b = AngularResidual.from_mesh(2.0 * np.ones(_ang_shape(m)), m)
        s = a + b
        d = b - a
        assert isinstance(s, AngularResidual) and isinstance(d, AngularResidual)
        assert np.all(s.values == 3.0)
        assert np.all(d.values == 1.0)
        # Originals untouched (frozen).
        assert np.all(a.values == 1.0)
        assert np.all(b.values == 2.0)

    def test_scalar_mul_div_neg(self) -> None:
        m = _slab_mesh()
        a = AngularResidual.from_mesh(np.full(_ang_shape(m), 4.0), m)
        assert np.all((3.0 * a).values == 12.0)
        assert np.all((a * 3.0).values == 12.0)
        assert np.all((a / 2.0).values == 2.0)
        assert np.all((-a).values == -4.0)
        assert isinstance(-a, AngularResidual)

    def test_mesh_binding_check(self) -> None:
        ma, mb = _slab_mesh(), _slab_mesh()
        a = AngularResidual.from_mesh(np.ones(_ang_shape(ma)), ma)
        b = AngularResidual.from_mesh(np.ones(_ang_shape(mb)), mb)
        with pytest.raises(ValueError, match="distinct SNMesh"):
            _ = a + b

    def test_frozen(self) -> None:
        m = _slab_mesh()
        r = AngularResidual.from_mesh(np.zeros(_ang_shape(m)), m)
        with pytest.raises(FrozenInstanceError):
            r.values = np.ones(_ang_shape(m))  # type: ignore[misc]


# ════════════════════════════════════════════════════════════════════
# ScalarResidual — construction + within-class algebra
# ════════════════════════════════════════════════════════════════════


class TestScalarResidual:
    def test_inherits_field_and_scalar_family(self) -> None:
        m = _slab_mesh()
        r = ScalarResidual.from_mesh(np.zeros(_sca_shape(m)), m)
        assert isinstance(r, Field)
        assert r.space.name == "scalar_residual"
        assert r.mesh is m

    def test_from_mesh_shape_and_metadata(self) -> None:
        m = _slab_mesh()
        r = ScalarResidual.from_mesh(np.ones(_sca_shape(m)), m)
        assert r.values.shape == _sca_shape(m)
        assert (r.ng, r.nx, r.ny) == (m.ng, m.nx, m.ny)

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="ScalarResidual.*does not match"):
            ScalarResidual.from_mesh(np.zeros((m.ng + 1, *m.spatial_shape)), m)

    def test_within_class_add_sub_closed(self) -> None:
        m = _slab_mesh()
        a = ScalarResidual.from_mesh(np.full(_sca_shape(m), 5.0), m)
        b = ScalarResidual.from_mesh(np.full(_sca_shape(m), 2.0), m)
        s = a + b
        d = a - b
        assert isinstance(s, ScalarResidual) and isinstance(d, ScalarResidual)
        assert np.all(s.values == 7.0)
        assert np.all(d.values == 3.0)

    def test_scalar_mul_div_neg(self) -> None:
        m = _slab_mesh()
        a = ScalarResidual.from_mesh(np.full(_sca_shape(m), 4.0), m)
        assert np.all((2.5 * a).values == 10.0)
        assert np.all((a / 4.0).values == 1.0)
        assert isinstance(-a, ScalarResidual)

    def test_mesh_binding_check(self) -> None:
        ma, mb = _slab_mesh(), _slab_mesh()
        a = ScalarResidual.from_mesh(np.ones(_sca_shape(ma)), ma)
        b = ScalarResidual.from_mesh(np.ones(_sca_shape(mb)), mb)
        with pytest.raises(ValueError, match="distinct SNMesh"):
            _ = a + b


# ════════════════════════════════════════════════════════════════════
# Cross-class rejection — the load-bearing B.3 invariant.
#
# Field's Layer-1 class-identity gate rejects EVERY cross-class
# operation BEFORE any value/space comparison, so "same units" never
# grants permission to combine. The balance that combines a residual
# with a source goes through a named composition (B.5), never a bare
# cross-class dunder.
# ════════════════════════════════════════════════════════════════════


class TestCrossClassRejection:
    def test_angular_residual_minus_angular_source_raises_same_units(
        self,
    ) -> None:
        r"""``AngularResidual - AngularSourceSink`` RAISES even though both
        carry :math:`1/(\mathrm{cm^3 \cdot s \cdot sr \cdot eV})`.

        Same shape, same units, DIFFERENT class — the dimensional-sin
        gate. Same units give permission to add in linear algebra; they
        do not give meaning.
        """
        m = _slab_mesh()
        res = AngularResidual.from_mesh(np.ones(_ang_shape(m)), m)
        src = AngularSourceSink.from_mesh(np.ones(_ang_shape(m)), m)
        # Identical shape and (conceptually) units — but distinct class.
        assert res.values.shape == src.values.shape
        with pytest.raises(TypeError, match="same-class"):
            _ = res - src  # type: ignore[operator]
        with pytest.raises(TypeError, match="same-class"):
            _ = src - res  # type: ignore[operator]

    def test_angular_residual_plus_angular_flux_raises(self) -> None:
        """``AngularResidual + AngularFlux`` RAISES (same shape, but a
        residual is a rate density and a flux is not)."""
        m = _slab_mesh()
        res = AngularResidual.from_mesh(np.ones(_ang_shape(m)), m)
        psi = AngularFlux.from_mesh(np.ones(_ang_shape(m)), m)
        with pytest.raises(TypeError, match="same-class"):
            _ = res + psi  # type: ignore[operator]

    def test_scalar_residual_minus_scalar_source_raises_same_units(
        self,
    ) -> None:
        """``ScalarResidual - ScalarSourceSink`` RAISES (same shape + units,
        distinct class)."""
        m = _slab_mesh()
        res = ScalarResidual.from_mesh(np.ones(_sca_shape(m)), m)
        src = ScalarSourceSink.from_mesh(np.ones(_sca_shape(m)), m)
        assert res.values.shape == src.values.shape
        with pytest.raises(TypeError, match="same-class"):
            _ = res - src  # type: ignore[operator]

    def test_scalar_residual_plus_scalar_flux_raises(self) -> None:
        m = _slab_mesh()
        res = ScalarResidual.from_mesh(np.ones(_sca_shape(m)), m)
        phi = ScalarFlux.from_mesh(np.ones(_sca_shape(m)), m)
        with pytest.raises(TypeError, match="same-class"):
            _ = res + phi  # type: ignore[operator]

    def test_cross_family_residual_add_raises(self) -> None:
        """``AngularResidual + ScalarResidual`` RAISES (different class
        AND different shape — the class gate fires first)."""
        m = _slab_mesh()
        ang = AngularResidual.from_mesh(np.ones(_ang_shape(m)), m)
        sca = ScalarResidual.from_mesh(np.ones(_sca_shape(m)), m)
        with pytest.raises(TypeError, match="same-class"):
            _ = ang + sca  # type: ignore[operator]

    def test_inner_product_cross_class_raises(self) -> None:
        """The cross-class gate also guards :meth:`Field.inner_product`."""
        m = _slab_mesh()
        res = AngularResidual.from_mesh(np.ones(_ang_shape(m)), m)
        src = AngularSourceSink.from_mesh(np.ones(_ang_shape(m)), m)
        with pytest.raises(TypeError, match="same-class"):
            res.inner_product(src)  # type: ignore[arg-type]


# ════════════════════════════════════════════════════════════════════
# 2-D smoke — the residual leaves work on a 2-D mesh as well.
# ════════════════════════════════════════════════════════════════════


class Test2DResiduals:
    def test_construction_shapes_2d(self) -> None:
        m = _2d_mesh()
        ang = AngularResidual.from_mesh(np.zeros(_ang_shape(m)), m)
        sca = ScalarResidual.from_mesh(np.zeros(_sca_shape(m)), m)
        assert ang.values.shape == _ang_shape(m)
        assert sca.values.shape == _sca_shape(m)

    def test_within_class_add_2d(self) -> None:
        m = _2d_mesh()
        rng = np.random.default_rng(31)
        a = AngularResidual.from_mesh(rng.standard_normal(_ang_shape(m)), m)
        b = AngularResidual.from_mesh(rng.standard_normal(_ang_shape(m)), m)
        c = a + b
        np.testing.assert_array_equal(c.values, a.values + b.values)


# ════════════════════════════════════════════════════════════════════
# from_balance — the B.5.1 named composition.
#
# Each residual leaf's bespoke factory: form the residual ``lhs − rhs``
# from two SAME-CLASS operands (the two sides of a balance) whose units
# match the residual's own. The role transition (e.g. AngularSourceSink
# operands → AngularResidual result) is the typed counterpart of the
# dimensional-sin rewire (B.5.2). Three guards: same-class operands,
# sr-exact units, same space + mesh.
# ════════════════════════════════════════════════════════════════════


class TestFromBalance:
    # ── Angular (bulk) ───────────────────────────────────────────────
    def test_angular_type_and_own_space(self) -> None:
        """Result is an ``AngularResidual`` on its OWN ``"angular_residual"``
        space — NOT the operands' ``"angular_source_sink"`` space."""
        m = _slab_mesh()
        lhs = AngularSourceSink.from_mesh(np.full(_ang_shape(m), 3.0), m)
        rhs = AngularSourceSink.from_mesh(np.full(_ang_shape(m), 1.0), m)
        r = AngularResidual.from_balance(lhs=lhs, rhs=rhs)
        assert isinstance(r, AngularResidual)
        assert r.space.name == "angular_residual"
        assert r.mesh is m

    def test_angular_sign_is_lhs_minus_rhs(self) -> None:
        m = _slab_mesh()
        lhs = AngularSourceSink.from_mesh(np.full(_ang_shape(m), 3.0), m)
        rhs = AngularSourceSink.from_mesh(np.full(_ang_shape(m), 1.0), m)
        r = AngularResidual.from_balance(lhs=lhs, rhs=rhs)
        np.testing.assert_array_equal(r.values, lhs.values - rhs.values)
        assert np.all(r.values == 2.0)

    def test_angular_positional_and_keyword_equivalent(self) -> None:
        m = _slab_mesh()
        rng = np.random.default_rng(7)
        a = AngularSourceSink.from_mesh(rng.standard_normal(_ang_shape(m)), m)
        b = AngularSourceSink.from_mesh(rng.standard_normal(_ang_shape(m)), m)
        np.testing.assert_array_equal(
            AngularResidual.from_balance(a, b).values,
            AngularResidual.from_balance(lhs=a, rhs=b).values,
        )

    def test_angular_result_addable_to_from_mesh_residual(self) -> None:
        """Load-bearing space-identity check: a ``from_balance`` residual and
        a ``from_mesh`` residual share one space, so they are mutually
        additive. (Would FAIL if ``from_balance`` reused the operands'
        ``"angular_source_sink"`` space.)"""
        m = _slab_mesh()
        lhs = AngularSourceSink.from_mesh(np.full(_ang_shape(m), 3.0), m)
        rhs = AngularSourceSink.from_mesh(np.full(_ang_shape(m), 1.0), m)
        from_bal = AngularResidual.from_balance(lhs=lhs, rhs=rhs)
        from_msh = AngularResidual.from_mesh(np.full(_ang_shape(m), 5.0), m)
        total = from_bal + from_msh  # MUST NOT raise (same space identity)
        assert isinstance(total, AngularResidual)
        assert np.all(total.values == 7.0)

    # ── Guards ───────────────────────────────────────────────────────
    def test_cross_class_operands_raise(self) -> None:
        """Operands of different class (source vs flux) → TypeError."""
        m = _slab_mesh()
        src = AngularSourceSink.from_mesh(np.ones(_ang_shape(m)), m)
        psi = AngularFlux.from_mesh(np.ones(_ang_shape(m)), m)
        with pytest.raises(TypeError, match="same-class operands"):
            AngularResidual.from_balance(lhs=src, rhs=psi)  # type: ignore[arg-type]

    def test_wrong_units_operands_raise(self) -> None:
        """Operands are a valid same-class pair (``AngularFlux``), but their
        units ``1/(cm²·s·sr)`` differ from ``AngularResidual``'s
        ``1/(cm³·s·sr)`` → TypeError (the sr-exact dimensional guard)."""
        m = _slab_mesh()
        a = AngularFlux.from_mesh(np.ones(_ang_shape(m)), m)
        b = AngularFlux.from_mesh(np.ones(_ang_shape(m)), m)
        with pytest.raises(TypeError, match="must match the residual's units"):
            AngularResidual.from_balance(lhs=a, rhs=b)  # type: ignore[arg-type]

    def test_mesh_mismatch_raises(self) -> None:
        ma, mb = _slab_mesh(), _slab_mesh()
        a = AngularSourceSink.from_mesh(np.ones(_ang_shape(ma)), ma)
        b = AngularSourceSink.from_mesh(np.ones(_ang_shape(mb)), mb)
        with pytest.raises(ValueError, match="distinct SNMesh"):
            AngularResidual.from_balance(lhs=a, rhs=b)

    # ── Scalar (bulk) ────────────────────────────────────────────────
    def test_scalar_type_sign_space(self) -> None:
        m = _slab_mesh()
        lhs = ScalarSourceSink.from_mesh(np.full(_sca_shape(m), 5.0), m)
        rhs = ScalarSourceSink.from_mesh(np.full(_sca_shape(m), 2.0), m)
        r = ScalarResidual.from_balance(lhs=lhs, rhs=rhs)
        assert isinstance(r, ScalarResidual)
        assert r.space.name == "scalar_residual"
        assert np.all(r.values == 3.0)

    def test_scalar_2d(self) -> None:
        m = _2d_mesh()
        rng = np.random.default_rng(13)
        lhs = ScalarSourceSink.from_mesh(rng.standard_normal(_sca_shape(m)), m)
        rhs = ScalarSourceSink.from_mesh(rng.standard_normal(_sca_shape(m)), m)
        r = ScalarResidual.from_balance(lhs=lhs, rhs=rhs)
        np.testing.assert_array_equal(r.values, lhs.values - rhs.values)


# ════════════════════════════════════════════════════════════════════
# Boundary locus — the uniform from_mesh factory + from_balance (the
# capability is minted in B.5.1; the matvec wiring is deferred to
# B.5.2 / #208). All boundary leaves share the SAME ``mesh.trace`` space.
# ════════════════════════════════════════════════════════════════════


class TestBoundaryFromMeshAndBalance:
    def test_from_mesh_roundtrip_on_trace(self) -> None:
        """``BoundaryField.from_mesh`` (B.5.1) packs a flat buffer onto the
        shared ``mesh.trace`` — the uniform construction surface bulk leaves
        already had."""
        m = _slab_mesh()
        n = BoundaryFlux.zeros_on(m).values.size
        vals = np.arange(float(n))
        bf = BoundaryFlux.from_mesh(vals, m)
        assert isinstance(bf, BoundaryFlux)
        assert bf.space is m.trace
        np.testing.assert_array_equal(bf.values, vals)

    def test_boundary_from_balance(self) -> None:
        m = _slab_mesh()
        n = BoundaryFlux.zeros_on(m).values.size
        lhs = BoundaryFlux.from_mesh(np.full(n, 4.0), m)
        rhs = BoundaryFlux.from_mesh(np.full(n, 1.0), m)
        r = BoundaryResidual.from_balance(lhs=lhs, rhs=rhs)
        assert isinstance(r, BoundaryResidual)
        assert r.space is m.trace  # shared trace space (all boundary leaves)
        assert np.all(r.values == 3.0)
