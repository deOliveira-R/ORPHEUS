r"""Foundation tests for :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`.

The typed angular spherical-harmonic frame's contract:

* :meth:`from_galerkin` upgrades a generic ``quadrature.angular_frame(L)``
  :class:`~orpheus.numerics.frame.GalerkinFrame` to the typed wrapper, sharing
  its basis + measure (table bit-identical);
* the role-polymorphic verbs :meth:`analyse` (:math:`M`) and :meth:`reconstruct`
  (:math:`R`) dispatch on the input carrier type, preserving the role
  (flux ↔ flux, source ↔ source) and changing the axis (angular ↔ moment);
* each verb is **bit-identical** to the generic ``np.ndarray`` face it wraps
  (``frame.analysis.apply`` / ``frame.reconstruction.apply``) — the typed seam
  adds carriers, never changes a number;
* a wrong-role carrier raises ``TypeError``.

Foundation mark: these verify a software invariant (the typed wrapper equals the
raw-face path) rather than a physics claim.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux
from orpheus.transport.frames import HarmonicFrame
from orpheus.transport.source_sinks import (
    AngularSourceSink,
    HarmonicMomentSourceSink,
)

from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

_L = 2


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
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
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _frame(m: SNMesh, L: int = _L) -> HarmonicFrame:
    return HarmonicFrame.from_galerkin(m.quad.angular_frame(L))


def _angular_values(m: SNMesh, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.standard_normal((m.quad.N, m.ng, *m.spatial_shape))


def _moment_values(m: SNMesh, seed: int, L: int = _L) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.standard_normal((L + 1, 2 * L + 1, m.ng, *m.spatial_shape))


# ── from_galerkin: the upgrade is a behaviour-identical specialization ──


class TestFromGalerkin:
    def test_is_a_galerkin_frame_sharing_basis_measure(self) -> None:
        m = _slab_mesh()
        gf = m.quad.angular_frame(_L)
        frame = HarmonicFrame.from_galerkin(gf)
        assert isinstance(frame, type(gf))  # Liskov: HarmonicFrame IS-A GalerkinFrame
        assert frame.basis is gf.basis
        assert frame.measure is gf.measure
        # the trial table (Yₗᵐ tabulation) is bit-identical
        np.testing.assert_array_equal(frame.table, gf.table)

    def test_inherited_faces_unchanged(self) -> None:
        """The generic ndarray faces are inherited untouched (0-ULP-safe)."""
        m = _slab_mesh()
        gf = m.quad.angular_frame(_L)
        frame = HarmonicFrame.from_galerkin(gf)
        vals = _angular_values(m, 1)
        np.testing.assert_array_equal(
            frame.analysis.apply(vals), gf.analysis.apply(vals),
        )


# ── analyse (M): role-preserving, angular → moment, bit-identical ──────


class TestAnalyse:
    def test_angular_flux_to_moment_flux_bit_identical(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        vals = _angular_values(m, 2)
        psi = AngularFlux.from_mesh(vals, m)
        moments = frame.analyse(psi)
        assert isinstance(moments, HarmonicMomentFlux)
        assert moments.mesh is m
        assert moments.L == _L
        np.testing.assert_array_equal(
            moments.values, frame.analysis.apply(vals),
        )

    def test_angular_source_to_moment_source_bit_identical(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        vals = _angular_values(m, 3)
        q = AngularSourceSink.from_mesh(vals, m)
        moments = frame.analyse(q)
        assert isinstance(moments, HarmonicMomentSourceSink)
        assert moments.mesh is m
        np.testing.assert_array_equal(
            moments.values, frame.analysis.apply(vals),
        )

    def test_wrong_carrier_raises(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        moment = HarmonicMomentFlux.from_mesh_and_L(_moment_values(m, 4), m, _L)
        with pytest.raises(TypeError, match="unsupported carrier"):
            frame.analyse(moment)


# ── reconstruct (R): role-preserving, moment → angular, bit-identical ──


class TestReconstruct:
    def test_moment_flux_to_angular_flux_bit_identical(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        vals = _moment_values(m, 5)
        phi = HarmonicMomentFlux.from_mesh_and_L(vals, m, _L)
        psi = frame.reconstruct(phi)
        assert isinstance(psi, AngularFlux)
        assert psi.mesh is m
        np.testing.assert_array_equal(
            psi.values, frame.reconstruction.apply(vals),
        )

    def test_moment_source_to_angular_source_bit_identical(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        vals = _moment_values(m, 6)
        q = HarmonicMomentSourceSink.from_mesh_and_L(vals, m, _L)
        out = frame.reconstruct(q)
        assert isinstance(out, AngularSourceSink)
        assert out.mesh is m
        np.testing.assert_array_equal(
            out.values, frame.reconstruction.apply(vals),
        )

    def test_wrong_carrier_raises(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        psi = AngularFlux.from_mesh(_angular_values(m, 7), m)
        with pytest.raises(TypeError, match="unsupported carrier"):
            frame.reconstruct(psi)


# ── role symmetry + 2-D smoke ──────────────────────────────────────────


class TestRolePolymorphism:
    def test_analyse_then_reconstruct_preserves_role_flux(self) -> None:
        """``reconstruct(analyse(flux))`` stays in the flux role (the band-
        limited projection Π = R∘M; not identity, but role-preserving)."""
        m = _slab_mesh()
        frame = _frame(m)
        psi = AngularFlux.from_mesh(_angular_values(m, 8), m)
        back = frame.reconstruct(frame.analyse(psi))
        assert isinstance(back, AngularFlux)
        # bit-identical to the raw composed faces
        expected = frame.reconstruction.apply(frame.analysis.apply(psi.values))
        np.testing.assert_array_equal(back.values, expected)

    def test_analyse_then_reconstruct_preserves_role_source(self) -> None:
        m = _slab_mesh()
        frame = _frame(m)
        q = AngularSourceSink.from_mesh(_angular_values(m, 9), m)
        back = frame.reconstruct(frame.analyse(q))
        assert isinstance(back, AngularSourceSink)

    def test_2d_analyse_reconstruct(self) -> None:
        m = _2d_mesh()
        frame = _frame(m)
        psi = AngularFlux.from_mesh(_angular_values(m, 10), m)
        moments = frame.analyse(psi)
        assert isinstance(moments, HarmonicMomentFlux)
        assert moments.values.shape == (_L + 1, 2 * _L + 1, m.ng, *m.spatial_shape)
        assert isinstance(frame.reconstruct(moments), AngularFlux)
