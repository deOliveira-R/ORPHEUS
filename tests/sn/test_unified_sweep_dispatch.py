"""Foundation tests for the unified sweep dispatch logic.

Round 2 of Wave D of the SN reshape campaign (Issue #161). These
tests pin the **dispatch contract** of :func:`orpheus.sn.sweep.transport_sweep`:

* It dispatches via
  :attr:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.requires_upstream_angular_state`
  (Wave B/D primitive), NOT by string-comparing
  ``sn_mesh.curvature``.  Slab + 2-D Cartesian → the Cartesian sweep;
  spherical + cylindrical → the curvilinear sweep.
* The 1-D cumprod fast path activates **only** when all three
  preconditions hold simultaneously: ``cell_update is
  DiamondDifference``, GL1D quadrature, isotropic source.  Any
  precondition failure routes through the 2-D wavefront sweep
  (handles 1-D as a special case of 2-D).

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
    _sweep_1d_cartesian,
    _sweep_1d_curvilinear,
    _sweep_2d_wavefront,
    transport_sweep,
)


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
    return SNMesh(mesh, quad)


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
    return SNMesh(mesh, quad)


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
    return SNMesh(mesh, quad)


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
    return SNMesh(mesh, quad)


# ═══════════════════════════════════════════════════════════════════════
# TestDispatchByReducedProperty
# ═══════════════════════════════════════════════════════════════════════

class TestDispatchByReducedProperty:
    """``transport_sweep`` dispatches via
    ``L.geometry.reduced.requires_upstream_angular_state``.

    The pre-Wave-D dispatch did string equality on ``sn_mesh.curvature``
    (``"spherical"``, ``"cylindrical"``, or ``None``).  Wave D Round 1
    exposed the ``ReducedStreamingOperator`` on
    :class:`SNMesh`; Wave D Round 2 (this issue) drives dispatch from
    the operator's ``requires_upstream_angular_state`` boolean — the
    geometry primitive that already encodes "does the sweep need an
    upstream half-angle ψ to apply redistribution?".
    """

    @pytest.mark.foundation
    def test_slab_routes_to_sweep_1d_cartesian(self, monkeypatch):
        """Slab (``requires_upstream_angular_state == False``) routes
        through ``_sweep_1d_cartesian``.
        """
        sn_mesh = _slab_sn_mesh()
        assert sn_mesh.reduced.requires_upstream_angular_state is False

        called = {"cartesian": 0, "curvilinear": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_cartesian(*args, **kwargs):
            called["cartesian"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_curvilinear(*args, **kwargs):
            called["curvilinear"] += 1
            return np.zeros((1,)), np.zeros((1,))

        monkeypatch.setattr(sweep_module, "_sweep_1d_cartesian", fake_cartesian)
        monkeypatch.setattr(sweep_module, "_sweep_1d_curvilinear", fake_curvilinear)

        nx, ng = 8, 1
        Q = np.zeros((nx, 1, ng))
        sig_t = np.ones((nx, 1, ng))
        transport_sweep(Q, sig_t, sn_mesh, {})

        assert called["cartesian"] == 1
        assert called["curvilinear"] == 0

    @pytest.mark.foundation
    def test_spherical_routes_to_sweep_1d_curvilinear(self, monkeypatch):
        """Spherical (``requires_upstream_angular_state == True``)
        routes through ``_sweep_1d_curvilinear``.
        """
        sn_mesh = _spherical_sn_mesh()
        assert sn_mesh.reduced.requires_upstream_angular_state is True

        called = {"cartesian": 0, "curvilinear": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_cartesian(*args, **kwargs):
            called["cartesian"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_curvilinear(*args, **kwargs):
            called["curvilinear"] += 1
            return np.zeros((1,)), np.zeros((1,))

        monkeypatch.setattr(sweep_module, "_sweep_1d_cartesian", fake_cartesian)
        monkeypatch.setattr(sweep_module, "_sweep_1d_curvilinear", fake_curvilinear)

        nx, ng = 8, 1
        Q = np.zeros((nx, 1, ng))
        sig_t = np.ones((nx, 1, ng))
        transport_sweep(Q, sig_t, sn_mesh, {})

        assert called["cartesian"] == 0
        assert called["curvilinear"] == 1

    @pytest.mark.foundation
    def test_cylindrical_routes_to_sweep_1d_curvilinear(self, monkeypatch):
        """Cylindrical (``requires_upstream_angular_state == True``)
        routes through ``_sweep_1d_curvilinear``.
        """
        sn_mesh = _cylindrical_sn_mesh()
        assert sn_mesh.reduced.requires_upstream_angular_state is True

        called = {"cartesian": 0, "curvilinear": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_cartesian(*args, **kwargs):
            called["cartesian"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_curvilinear(*args, **kwargs):
            called["curvilinear"] += 1
            return np.zeros((1,)), np.zeros((1,))

        monkeypatch.setattr(sweep_module, "_sweep_1d_cartesian", fake_cartesian)
        monkeypatch.setattr(sweep_module, "_sweep_1d_curvilinear", fake_curvilinear)

        nx, ng = 8, 1
        Q = np.zeros((nx, 1, ng))
        sig_t = np.ones((nx, 1, ng))
        transport_sweep(Q, sig_t, sn_mesh, {})

        assert called["cartesian"] == 0
        assert called["curvilinear"] == 1

    @pytest.mark.foundation
    def test_2d_cartesian_routes_to_sweep_2d_wavefront(self, monkeypatch):
        """2-D Cartesian (``sn_mesh.reduced is None``, ``ny > 1``) routes
        through ``_sweep_2d_wavefront``.

        Issue #196 Step 2.5: the top-level dispatch now separates 1-D
        from 2-D Cartesian directly (``sn_mesh.ny == 1`` →
        :func:`_sweep_1d_cartesian`; ``ny > 1`` →
        :func:`_sweep_2d_wavefront`).  Pre-Step-2.5 these were both
        called through a ``_cartesian_sweep`` indirection that did
        the 1D-vs-2D split internally.
        """
        sn_mesh = _2d_sn_mesh()
        assert sn_mesh.reduced is None
        assert sn_mesh.ny > 1

        called = {"cartesian": 0, "wavefront": 0, "curvilinear": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_cartesian(*args, **kwargs):
            called["cartesian"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            return np.zeros((1,)), np.zeros((1,))

        def fake_curvilinear(*args, **kwargs):
            called["curvilinear"] += 1
            return np.zeros((1,)), np.zeros((1,))

        monkeypatch.setattr(sweep_module, "_sweep_1d_cartesian", fake_cartesian)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)
        monkeypatch.setattr(sweep_module, "_sweep_1d_curvilinear", fake_curvilinear)

        nx, ny, ng = 4, 4, 1
        Q = np.zeros((nx, ny, ng))
        sig_t = np.ones((nx, ny, ng))
        transport_sweep(Q, sig_t, sn_mesh, {})

        assert called["cartesian"] == 0
        assert called["wavefront"] == 1
        assert called["curvilinear"] == 0


# ═══════════════════════════════════════════════════════════════════════
# TestCartesianDispatch1Dvs2D — Issue #196 Step 2.5
# ═══════════════════════════════════════════════════════════════════════

class TestCartesianDispatch1Dvs2D:
    """Cartesian dispatch separates 1-D from 2-D (Issue #196 Step 2.5).

    Pre-Step-2.5 the dispatch had a 1-D cumprod fast path gated by
    three preconditions (DD strategy, GL1D quadrature, isotropic
    source).  Step 2.5 retired the cumprod path and the gating logic
    along with it: 1-D Cartesian (``sn_mesh.ny == 1``) routes
    unconditionally through :func:`_sweep_1d_cartesian` (the unified
    fold-over-DAG-visits skeleton); 2-D Cartesian routes through
    :func:`_sweep_2d_wavefront`.  Cell-update strategy, quadrature
    shape, and source anisotropy are all handled INSIDE the chosen
    sweep — they no longer affect dispatch.
    """

    @pytest.mark.foundation
    def test_1d_cartesian_routes_to_sweep_1d_cartesian(self, monkeypatch):
        """1-D Cartesian (``sn_mesh.ny == 1``) → ``_sweep_1d_cartesian``."""
        sn_mesh = _slab_sn_mesh()
        assert sn_mesh.ny == 1
        assert sn_mesh.reduced.requires_upstream_angular_state is False

        called = {"cartesian": 0, "wavefront": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_cartesian(*args, **kwargs):
            called["cartesian"] += 1
            quad_N = sn_mesh.quad.N
            return (np.zeros((quad_N, sn_mesh.nx, 1, 1)),
                    np.zeros((sn_mesh.nx, 1, 1)))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            return (np.zeros((1,)), np.zeros((1,)))

        monkeypatch.setattr(sweep_module, "_sweep_1d_cartesian", fake_cartesian)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)

        nx, ng = 8, 1
        Q = np.zeros((nx, 1, ng))
        sig_t = np.ones((nx, 1, ng))
        transport_sweep(Q, sig_t, sn_mesh, {})

        assert called["cartesian"] == 1
        assert called["wavefront"] == 0

    @pytest.mark.foundation
    def test_1d_cartesian_anisotropic_routes_to_sweep_1d_cartesian(
        self, monkeypatch,
    ):
        """1-D + anisotropic source still routes to ``_sweep_1d_cartesian``.

        Step 2.5 retired the fast-path gating: the unified
        ``_sweep_1d_cartesian`` body handles ``Q_aniso`` natively
        (per-ordinate anisotropic source folded into ``QV`` inside
        the per-ordinate loop).
        """
        sn_mesh = _slab_sn_mesh()
        assert sn_mesh.ny == 1

        called = {"cartesian": 0, "wavefront": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_cartesian(*args, **kwargs):
            called["cartesian"] += 1
            quad_N = sn_mesh.quad.N
            return (np.zeros((quad_N, sn_mesh.nx, 1, 1)),
                    np.zeros((sn_mesh.nx, 1, 1)))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            return (np.zeros((1,)), np.zeros((1,)))

        monkeypatch.setattr(sweep_module, "_sweep_1d_cartesian", fake_cartesian)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)

        nx, ng = 8, 1
        Q = np.zeros((nx, 1, ng))
        sig_t = np.ones((nx, 1, ng))
        Q_aniso = np.zeros((sn_mesh.quad.N, nx, 1, ng))  # NOT None
        transport_sweep(Q, sig_t, sn_mesh, {}, Q_aniso=Q_aniso)

        assert called["cartesian"] == 1
        assert called["wavefront"] == 0

    @pytest.mark.foundation
    def test_2d_cartesian_routes_to_sweep_2d_wavefront(self, monkeypatch):
        """2-D Cartesian (``sn_mesh.ny > 1``) → ``_sweep_2d_wavefront``."""
        sn_mesh = _2d_sn_mesh()
        assert sn_mesh.ny > 1
        assert sn_mesh.reduced is None

        called = {"cartesian": 0, "wavefront": 0}
        from orpheus.sn import sweep as sweep_module

        def fake_cartesian(*args, **kwargs):
            called["cartesian"] += 1
            return (np.zeros((1,)), np.zeros((1,)))

        def fake_wavefront(*args, **kwargs):
            called["wavefront"] += 1
            quad_N = sn_mesh.quad.N
            return (np.zeros((quad_N, sn_mesh.nx, sn_mesh.ny, 1)),
                    np.zeros((sn_mesh.nx, sn_mesh.ny, 1)))

        monkeypatch.setattr(sweep_module, "_sweep_1d_cartesian", fake_cartesian)
        monkeypatch.setattr(sweep_module, "_sweep_2d_wavefront", fake_wavefront)

        nx, ny, ng = 4, 4, 1
        Q = np.zeros((nx, ny, ng))
        sig_t = np.ones((nx, ny, ng))
        transport_sweep(Q, sig_t, sn_mesh, {})

        assert called["cartesian"] == 0
        assert called["wavefront"] == 1


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
        sn_mesh = SNMesh(mesh, quad, cell_update=custom)
        assert sn_mesh.cell_update is custom
