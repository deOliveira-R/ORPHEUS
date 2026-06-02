r"""O.1 (Issue #208 / Wave O) — SN operator-leaf block-role tags.

Pins that every SN operator leaf advertises the correct block-role
classification on the 2×2 bulk/boundary block structure:

* Collision ``C``, scattering ``S``, fission ``F`` → :attr:`BlockRole.BULK`
  (``A_bb`` only; read/write the bulk flux, no boundary action).
* Streaming ``L`` → :attr:`BlockRole.FULL` (couples bulk ↔ boundary —
  reads the inflow trace to seed the sweep, writes the outflow trace).
* ``(L + C)`` :class:`InvertibleOperator` → :attr:`BlockRole.FULL`
  (inherits the streaming operand's coupling).

The classification MECHANISM (enum, markers, partition, None-default) is
pinned in ``tests/numerics/test_operator_protocols.py``; here we pin the
SN tagging itself.

For ``C``, ``L`` and ``(L+C)`` we construct instances and assert
``isinstance`` + exclusivity. For ``S`` and ``F`` (which need a
``MaterialXSField``) we assert the class-level tag: because ``block_role``
is a class attribute and the marker metaclass reads the *instance's*
``block_role``, ``isinstance(s_instance, BulkOperator)`` is by
construction equivalent to ``ScatteringOperator.block_role is
BlockRole.BULK`` — the class-level assertion fully determines the
instance membership, so it is not a coverage gap.

The BOUNDARY role has no first-class leaf until O.4 (the boundary law is
absorbed in ``L`` today), so it is not exercised here.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.numerics.operator import BlockRole, BulkOperator, FullOperator
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.fission import FissionOperator
from orpheus.sn.operator import (
    CollisionOperator,
    InvertibleOperator,
    StreamingOperator,
)
from orpheus.sn.scattering import ScatteringOperator
from tests.sn._test_helpers import placeholder_materials

pytestmark = [pytest.mark.foundation]


def _slab_mesh(nx: int = 4, n_ord: int = 4, ng: int = 1) -> SNMesh:
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.vacuum, BC.vacuum),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=nx),))
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


class TestBulkLeaves:
    def test_collision_is_bulk(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        C = CollisionOperator(sn, sigma_t)
        assert C.block_role is BlockRole.BULK
        assert isinstance(C, BulkOperator)
        assert not isinstance(C, FullOperator)

    def test_scattering_is_bulk(self) -> None:
        assert ScatteringOperator.block_role is BlockRole.BULK
        assert ScatteringOperator.block_role is not BlockRole.FULL

    def test_fission_is_bulk(self) -> None:
        assert FissionOperator.block_role is BlockRole.BULK
        assert FissionOperator.block_role is not BlockRole.FULL


class TestFullLeaves:
    def test_streaming_is_full(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        L = StreamingOperator(sn, sigma_t)
        assert L.block_role is BlockRole.FULL
        assert isinstance(L, FullOperator)
        assert not isinstance(L, BulkOperator)

    def test_invertible_L_plus_C_is_full(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, sn.nx, sn.ny))
        composite = StreamingOperator(sn, sigma_t) + CollisionOperator(sn, sigma_t)
        assert isinstance(composite, InvertibleOperator)
        assert composite.block_role is BlockRole.FULL
        assert isinstance(composite, FullOperator)
        assert not isinstance(composite, BulkOperator)
