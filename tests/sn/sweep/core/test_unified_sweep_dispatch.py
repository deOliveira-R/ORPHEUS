"""Foundation tests for the sweep-strategy dispatch.

Round 2 of Wave D of the SN reshape campaign (Issue #161); migrated by
the **sweep-strategy carve** (C3.4/C3.5, plan ``sn_sweep_strategy.md``).

Historically ``transport_sweep`` chose the 1-D vs 2-D sweep body with a
scattered ``sn_mesh.reduced is not None`` branch, and these tests spied
on the chosen ``_sweep_*`` function being *called* (monkeypatching the
module-level name).  The carve replaced that branch with a first-class,
selectable :class:`~orpheus.sn.loss_representation.LossRepresentation`:
:func:`~orpheus.sn.loss_representation.default_for` picks the strategy, and
``transport_sweep`` delegates to ``strategy.sweep``.  The spy mechanism no
longer applies (the strategy holds its own reference to the sweep body),
so the dispatch contract is now pinned at its single source of truth:

* the SELECTION — ``default_for(mesh)`` returns the right strategy class:
  ALL 1-D meshes (slab, sphere, cylinder; ``is_1d``) →
  :class:`~orpheus.sn.loss_representation.CumprodScan`; 2-D Cartesian
  (``is_cartesian and ndim == 2``) →
  :class:`~orpheus.sn.loss_representation.MovingFrontierWindow`;
* the ROUTING — ``transport_sweep`` delegates to
  ``default_for(mesh).sweep(...)`` exactly once.

``-O``-safe (vv Mode 8): every gate is a ``np.testing`` / ``pytest.fail``
function call, NOT a bare ``assert`` (which the canonical ``python -O``
invocation would strip in any non-rewritten helper).

These are **software-contract** tests — the L1 transport math is
verified by the regression snapshots at
``tests/sn/regression/snapshots/`` (gate the bit-identical contract)
and the Wave C ``DiamondDifference`` hand-calc tests at
``tests/sn/spatial/test_diamond.py`` (gate the per-cell algebra).

Tagged ``@pytest.mark.foundation`` because:

* The dispatch is purely structural — no transport equation identity is
  being verified here.
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
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.loss_representation import transport_sweep
from orpheus.sn.loss_representation import (
    CumprodScan,
    MovingFrontierWindow,
    default_for,
)
from tests.sn._test_helpers import placeholder_materials
from orpheus.transport.source_sinks import AngularSourceSink
from orpheus.transport.fields.boundary_flux import BoundaryFlux


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
    quad = Quadrature.gauss_legendre(n_ordinates=8)
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
    quad = Quadrature.gauss_legendre(n_ordinates=8)
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
    quad = Quadrature.level_symmetric(sn_order=4)
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
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials())


# ═══════════════════════════════════════════════════════════════════════
# TestDispatchSelectsStrategy — the SELECTION half of the contract
# ═══════════════════════════════════════════════════════════════════════

class TestDispatchSelectsStrategy:
    """``default_for(mesh)`` selects the geometry-correct sweep strategy.

    The single source of truth that replaced the scattered
    ``reduced is not None`` branch: ALL 1-D meshes (slab, sphere, cylinder)
    → :class:`CumprodScan`; 2-D Cartesian → :class:`MovingFrontierWindow`
    (the production optimization; the full-field oracle is explicit-select
    only).  ``pytest.fail`` on mismatch fires under ``python -O``.
    """

    @pytest.mark.foundation
    def test_slab_selects_cumprod_scan(self):
        """Slab (1-D Cartesian) → CumprodScan."""
        sn_mesh = _slab_sn_mesh()
        if sn_mesh.reduced is None:
            pytest.fail("slab fixture unexpectedly has reduced is None")
        strategy = default_for(sn_mesh)
        if not isinstance(strategy, CumprodScan):
            pytest.fail(
                f"slab → {type(strategy).__name__}, expected CumprodScan"
            )

    @pytest.mark.foundation
    def test_spherical_selects_cumprod_scan(self):
        """Spherical (1-D curvilinear) → CumprodScan (same as slab)."""
        sn_mesh = _spherical_sn_mesh()
        if sn_mesh.reduced is None:
            pytest.fail("spherical fixture unexpectedly has reduced is None")
        strategy = default_for(sn_mesh)
        if not isinstance(strategy, CumprodScan):
            pytest.fail(
                f"sphere → {type(strategy).__name__}, expected CumprodScan"
            )

    @pytest.mark.foundation
    def test_cylindrical_selects_cumprod_scan(self):
        """Cylindrical (1-D curvilinear) → CumprodScan (same as slab/sphere)."""
        sn_mesh = _cylindrical_sn_mesh()
        if sn_mesh.reduced is None:
            pytest.fail("cylindrical fixture unexpectedly has reduced is None")
        strategy = default_for(sn_mesh)
        if not isinstance(strategy, CumprodScan):
            pytest.fail(
                f"cylinder → {type(strategy).__name__}, expected CumprodScan"
            )

    @pytest.mark.foundation
    def test_2d_cartesian_selects_moving_frontier_window(self):
        """2-D Cartesian → MovingFrontierWindow (the production optimization)."""
        sn_mesh = _2d_sn_mesh()
        if sn_mesh.reduced is not None:
            pytest.fail("2-D fixture unexpectedly has reduced is not None")
        strategy = default_for(sn_mesh)
        if not isinstance(strategy, MovingFrontierWindow):
            pytest.fail(
                f"2-D Cartesian → {type(strategy).__name__}, "
                f"expected MovingFrontierWindow"
            )


# ═══════════════════════════════════════════════════════════════════════
# TestTransportSweepDelegatesToStrategy — the ROUTING half of the contract
# ═══════════════════════════════════════════════════════════════════════

class TestTransportSweepDelegatesToStrategy:
    """``transport_sweep`` routes through ``default_for(mesh).sweep`` once.

    The ROUTING half of the dispatch contract (the SELECTION half is
    :class:`TestDispatchSelectsStrategy`).  ``transport_sweep`` does a lazy
    ``from .loss_representation import default_for``, so patching
    ``loss_representation.default_for`` is seen at call time — the spy confirms the
    dispatcher *delegates to the selected strategy* rather than re-deciding
    the branch itself.  Geometry-agnostic: the same delegation holds for the
    1-D scan and the 2-D window.
    """

    @pytest.mark.foundation
    @pytest.mark.parametrize(
        "mesh_factory", [_slab_sn_mesh, _2d_sn_mesh], ids=["slab", "cart2d"],
    )
    def test_delegates_to_selected_strategy(self, monkeypatch, mesh_factory):
        """transport_sweep calls the selected strategy's ``sweep`` exactly once."""
        sn_mesh = mesh_factory()
        import orpheus.sn.loss_representation as loss_representation

        selected = type(default_for(sn_mesh)).__name__
        calls = {"sweep": 0}
        N, ng, nx, ny = sn_mesh.quad.N, sn_mesh.ng, sn_mesh.nx, sn_mesh.ny

        class _SpyStrategy:
            def sweep(self, *args, **kwargs):
                calls["sweep"] += 1
                return (np.zeros((N, ng, nx, ny)), np.zeros((ng, nx, ny)))

        monkeypatch.setattr(
            loss_representation, "default_for", lambda mesh: _SpyStrategy(),
        )

        # 1-D meshes carry (ng, nx) Σ_t; 2-D Cartesian carries (ng, nx, ny).
        sig_t = (
            np.ones((ng, nx)) if sn_mesh.reduced is not None
            else np.ones((ng, nx, ny))
        )
        source = AngularSourceSink.zeros_on(sn_mesh)
        transport_sweep(source, sig_t, sn_mesh, BoundaryFlux.zeros_on(sn_mesh))

        if calls["sweep"] != 1:
            pytest.fail(
                f"transport_sweep delegated to strategy.sweep "
                f"{calls['sweep']} times (selected {selected}), expected "
                f"exactly 1"
            )


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
        if not isinstance(sn_mesh.cell_update, DiamondDifference):
            pytest.fail(
                f"default cell_update is {type(sn_mesh.cell_update).__name__}, "
                f"expected DiamondDifference"
            )

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
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        sn_mesh = SNMesh(mesh, quad, placeholder_materials(), cell_update=custom)
        if sn_mesh.cell_update is not custom:
            pytest.fail("explicit cell_update was not stored on the mesh")
