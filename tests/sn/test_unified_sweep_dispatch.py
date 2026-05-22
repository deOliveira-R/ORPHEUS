"""Foundation tests for the unified sweep dispatch logic.

Round 2 of Wave D of the SN reshape campaign (Issue #161).  Issue
#196 Phase G Step 2.5c collapses ``_sweep_1d_cartesian`` and
``_sweep_1d_curvilinear`` into ONE ``_sweep_1d_unified`` body driven
by the precomputed two-stratum cache (geometry is data on
:class:`GeometryCoefficients`, not control-flow at the dispatch
layer).  These tests now pin the simplified dispatch contract of
:func:`orpheus.sn.sweep.transport_sweep`:

* 1-D meshes (``sn_mesh.reduced is not None``) — slab, sphere,
  cylinder — ALL route to :func:`_sweep_1d_unified`.
* 2-D Cartesian (``sn_mesh.reduced is None``) routes to
  :func:`_sweep_2d_wavefront` (anti-diagonal scheduling; Step 2.6 Q2
  defers the 2D wavefront unification — see the plan).

These are **software-contract** tests — the L1 transport math is
verified by the regression snapshots at
``tests/sn/regression/snapshots/`` (gate the bit-identical contract)
and the Wave C ``DiamondDifference`` hand-calc tests at
``tests/sn/spatial/test_diamond.py`` (gate the per-cell algebra).

Tagged ``@pytest.mark.foundation`` because:

* The dispatch math is purely structural — no transport equation
  identity is being verified here.
* L1 transport accuracy is verified transitively via the existing
  MMS suite at ``tests/sn/l1_analytical/``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    Mesh2D,
)
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D, LevelSymmetricSN
from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.sweep import (
    _sweep_1d_unified,
    _sweep_2d_wavefront,
    transport_sweep,
)
from tests.sn._test_helpers import placeholder_materials


# ═══════════════════════════════════════════════════════════════════════
# Mesh fixtures
# ═══════════════════════════════════════════════════════════════════════

def _slab_sn_mesh(nx: int = 8, length: float = 1.0) -> SNMesh:
    """Slab SNMesh with vacuum BCs and Gauss-Legendre 1D quadrature."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=8)
    return SNMesh(mesh, quad, placeholder_materials())


def _spherical_sn_mesh(nx: int = 8, radius: float = 1.0) -> SNMesh:
    """Spherical SNMesh with reflective inner / vacuum outer BCs."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = GaussLegendre1D.create(n_ordinates=8)
    return SNMesh(mesh, quad, placeholder_materials())


def _cylindrical_sn_mesh(nx: int = 8, radius: float = 1.0) -> SNMesh:
    """Cylindrical SNMesh with reflective inner / vacuum outer BCs."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, radius, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials())


def _2d_sn_mesh(nx: int = 4, ny: int = 4) -> SNMesh:
    """2-D Cartesian SNMesh with vacuum BCs and Level-Symmetric S4."""
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"),
        bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"),
        bc_ymax=BC("vacuum"),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials())


# ═══════════════════════════════════════════════════════════════════════
# TestDispatchByReducedProperty
# ═══════════════════════════════════════════════════════════════════════

class TestDispatchByReducedProperty:
    """``transport_sweep`` dispatches by ``sn_mesh.reduced is not None``.

    Issue #196 Phase G Step 2.5c collapses the historical
    ``_sweep_1d_cartesian`` and ``_sweep_1d_curvilinear`` into ONE
    ``_sweep_1d_unified`` body (the two-stratum cache abstracts the
    geometry difference).  Dispatch is therefore: ALL 1-D meshes (slab,
    sphere, cylinder; ``sn_mesh.reduced is not None``) route to
    ``_sweep_1d_unified``; 2-D Cartesian (``sn_mesh.reduced is None``)
    routes to ``_sweep_2d_wavefront``.
    """

    @pytest.mark.foundation
    def test_slab_routes_to_sweep_1d_unified(self, monkeypatch):
        """Slab routes through ``_sweep_1d_unified``."""
        sn_mesh = _slab_sn_mesh()
        assert sn_mesh.reduced is not None

        called = {"unified": 0, "wavefront": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_unified(*args, **kwargs):
            called["unified"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            return np.zeros((1,)), np.zeros((1,))

        monkeypatch.setattr(sweep_module, "_sweep_1d_unified", fake_unified)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)

        # PR-INDEX-5 principled (ng, nx, ny); PR-TYPED-4 typed source.
        ng = sn_mesh.ng
        Q = sn_mesh.zeros_per_ordinate_source()
        sig_t = np.ones((ng, sn_mesh.nx, 1))
        transport_sweep(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())

        assert called["unified"] == 1
        assert called["wavefront"] == 0

    @pytest.mark.foundation
    def test_spherical_routes_to_sweep_1d_unified(self, monkeypatch):
        """Spherical routes through ``_sweep_1d_unified`` (same as slab)."""
        sn_mesh = _spherical_sn_mesh()
        assert sn_mesh.reduced is not None

        called = {"unified": 0, "wavefront": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_unified(*args, **kwargs):
            called["unified"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            return np.zeros((1,)), np.zeros((1,))

        monkeypatch.setattr(sweep_module, "_sweep_1d_unified", fake_unified)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)

        ng = sn_mesh.ng
        Q = sn_mesh.zeros_per_ordinate_source()
        sig_t = np.ones((ng, sn_mesh.nx, 1))
        transport_sweep(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())

        assert called["unified"] == 1
        assert called["wavefront"] == 0

    @pytest.mark.foundation
    def test_cylindrical_routes_to_sweep_1d_unified(self, monkeypatch):
        """Cylindrical routes through ``_sweep_1d_unified`` (same as slab/sphere)."""
        sn_mesh = _cylindrical_sn_mesh()
        assert sn_mesh.reduced is not None

        called = {"unified": 0, "wavefront": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_unified(*args, **kwargs):
            called["unified"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            return np.zeros((1,)), np.zeros((1,))

        monkeypatch.setattr(sweep_module, "_sweep_1d_unified", fake_unified)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)

        ng = sn_mesh.ng
        Q = sn_mesh.zeros_per_ordinate_source()
        sig_t = np.ones((ng, sn_mesh.nx, 1))
        transport_sweep(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())

        assert called["unified"] == 1
        assert called["wavefront"] == 0

    @pytest.mark.foundation
    def test_2d_cartesian_routes_to_sweep_2d_wavefront(self, monkeypatch):
        """2-D Cartesian (``sn_mesh.reduced is None``) routes to ``_sweep_2d_wavefront``."""
        sn_mesh = _2d_sn_mesh()
        assert sn_mesh.reduced is None

        called = {"unified": 0, "wavefront": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_unified(*args, **kwargs):
            called["unified"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            return np.zeros((1,)), np.zeros((1,))

        monkeypatch.setattr(sweep_module, "_sweep_1d_unified", fake_unified)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)

        ng = sn_mesh.ng
        Q = sn_mesh.zeros_per_ordinate_source()
        sig_t = np.ones((ng, sn_mesh.nx, sn_mesh.ny))
        transport_sweep(Q, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())

        assert called["unified"] == 0
        assert called["wavefront"] == 1


# ═══════════════════════════════════════════════════════════════════════
# TestUnifiedDispatch1Dvs2D — Issue #196 Step 2.5c
# ═══════════════════════════════════════════════════════════════════════

class TestUnifiedDispatch1Dvs2D:
    """1-D / 2-D dispatch contract under the Step 2.5c unification.

    Step 2.5c collapses both Cartesian and curvilinear 1-D paths into
    one body.  The dispatch contract is: ALL 1-D meshes go to
    ``_sweep_1d_unified``; 2-D meshes go to ``_sweep_2d_wavefront``.
    Cell-update strategy, quadrature shape, and source anisotropy are
    all handled INSIDE the chosen sweep — they no longer affect
    dispatch.
    """

    @pytest.mark.foundation
    def test_1d_unified_handles_anisotropic_source(self, monkeypatch):
        """1-D + anisotropic source still routes to ``_sweep_1d_unified``."""
        sn_mesh = _slab_sn_mesh()
        assert sn_mesh.ny == 1

        called = {"unified": 0, "wavefront": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_unified(*args, **kwargs):
            called["unified"] += 1
            quad_N = sn_mesh.quad.N
            return (np.zeros((quad_N, sn_mesh.nx, 1, 1)),
                    np.zeros((sn_mesh.nx, 1, 1)))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            return (np.zeros((1,)), np.zeros((1,)))

        monkeypatch.setattr(sweep_module, "_sweep_1d_unified", fake_unified)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)

        # R-1 Step 4 A1 — single per-ordinate source.
        ng = sn_mesh.ng
        source = sn_mesh.zeros_per_ordinate_source()
        sig_t = np.ones((ng, sn_mesh.nx, 1))
        transport_sweep(source, sig_t, sn_mesh, sn_mesh.zeros_boundary_flux())

        assert called["unified"] == 1
        assert called["wavefront"] == 0


# ═══════════════════════════════════════════════════════════════════════
# TestDefaultCellUpdate
# ═══════════════════════════════════════════════════════════════════════

class TestDefaultCellUpdate:
    """``SNMesh.cell_update`` defaults to :class:`DiamondDifference`.

    The default is what guarantees bit-identity with the pre-Wave-D
    sweep — DD's per-cell math is a bit-identical extraction of the
    inlined sweep's algebra.  Wave C-extension will let users pass
    LD / EC / Step strategies via the constructor argument.
    """

    @pytest.mark.foundation
    def test_default_is_diamond_difference(self):
        """No ``cell_update`` argument → defaults to DD."""
        sn_mesh = _slab_sn_mesh()
        assert isinstance(sn_mesh.cell_update, DiamondDifference)

    @pytest.mark.foundation
    def test_explicit_cell_update_honored(self):
        """User-passed ``cell_update`` is stored on the mesh."""
        custom = DiamondDifference()  # the only strategy that ships in Wave D
        mesh = Mesh1D(
            edges=np.linspace(0.0, 1.0, 9),
            mat_ids=np.zeros(8, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"),
            bc_right=BC("vacuum"),
        )
        quad = GaussLegendre1D.create(n_ordinates=8)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials(), cell_update=custom)
        assert sn_mesh.cell_update is custom
