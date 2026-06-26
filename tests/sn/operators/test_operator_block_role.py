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

The BOUNDARY role tags the realized boundary laws
(:meth:`SNBoundaryRealizer.realize` outputs — vacuum / reflective /
white / albedo / periodic carry :attr:`BlockRole.BOUNDARY`; the rank-0
affine ``PrescribedInflow`` source does NOT — it is the boundary
*source* ``q.boundary``, not a linear ``B``). Pinned in
``TestBoundaryLeaves`` below, both on the raw realizer output and on the
mesh-wired ``sn_mesh.bc`` entries (the ``_BoundBoundaryOperator`` shim forwards
the inner op's role). The O.4a.1-γ tagging lands the role; O.4a.2 wires
``B`` into ``(L_full + C − S − F − B)`` as a sibling of ``L``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    NoSource,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.numerics.operator import (
    BlockRole,
    BoundaryOperator,
    BulkOperator,
    FullOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace
from orpheus.sn.geometry import SNMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.sn.boundary_operator import SNBoundaryOperator
from orpheus.sn.operator import (
    InvertibleOperator,
    StreamingOperator,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.operators.scattering import ScatteringOperator
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
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        C = MultiplicationOperator.from_mesh(sigma_t, sn)
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
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        L = StreamingOperator(sn)
        assert L.block_role is BlockRole.FULL
        assert isinstance(L, FullOperator)
        assert not isinstance(L, BulkOperator)

    def test_invertible_L_plus_C_is_full(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        composite = StreamingOperator(sn) + MultiplicationOperator.from_mesh(sigma_t, sn)
        assert isinstance(composite, InvertibleOperator)
        assert composite.block_role is BlockRole.FULL
        assert isinstance(composite, FullOperator)
        assert not isinstance(composite, BulkOperator)


# ─────────────────────────────────────────────────────────────────────
# Boundary leaves (Wave O step O.4a.1-γ)
# ─────────────────────────────────────────────────────────────────────


def _boundary_method_space(n_ord: int = 4) -> SNMethodSpace:
    """A method space carrying ``inflow_indices`` (needed by the vacuum
    branch; harmless for the others — they read only ``quadrature``)."""
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    inflow = np.flatnonzero(quad.mu_x > 0)
    return SNMethodSpace(quadrature=quad, face="xmin", inflow_indices=inflow)


# Every LINEAR boundary law + a representative amplitude per kind. Each
# realizes to a BOUNDARY-block operator (``A_ss`` only). albedo=0/1 hit the
# Zero/Identity fast paths; α∉{0,1} hits the ScaledOperator path — all three
# must carry the role.
_LINEAR_LAWS = {
    "vacuum": VacuumInflow(),
    "reflective_albedo1": ReflectiveBoundary(axis="x", albedo=1.0),
    "reflective_albedo07": ReflectiveBoundary(axis="x", albedo=0.7),
    "white_albedo1": WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0),
    "white_albedo05": WhiteBoundary(axis="x", outward_sign=+1, albedo=0.5),
    "albedo_0": AlbedoBoundary(albedo=0.0),
    "albedo_1": AlbedoBoundary(albedo=1.0),
    "albedo_03": AlbedoBoundary(albedo=0.3),
    "periodic": PeriodicBoundary(),
}


class TestBoundaryLeaves:
    """The realized boundary laws advertise :attr:`BlockRole.BOUNDARY`."""

    @pytest.mark.parametrize("law_id", list(_LINEAR_LAWS))
    def test_realized_bc_advertises_boundary_operator(self, law_id: str) -> None:
        law = _LINEAR_LAWS[law_id]
        op = SNBoundaryRealizer().realize(law, _boundary_method_space())
        assert op.block_role is BlockRole.BOUNDARY
        assert isinstance(op, BoundaryOperator)
        # Exclusivity (vv L11): a BOUNDARY leaf is neither bulk nor full.
        assert not isinstance(op, BulkOperator)
        assert not isinstance(op, FullOperator)

    def test_prescribed_inflow_is_not_boundary_operator(self) -> None:
        """The rank-0 affine ``PrescribedInflow`` source is ``q.boundary``,
        NOT a linear boundary operator ``B`` — it carries no block role."""
        op = SNBoundaryRealizer().realize(
            PrescribedInflow(source=NoSource()), _boundary_method_space(),
        )
        assert getattr(op, "block_role", None) is None
        assert not isinstance(op, BoundaryOperator)
        assert not isinstance(op, BulkOperator)
        assert not isinstance(op, FullOperator)

    def test_mesh_bc_forwards_boundary_role(self) -> None:
        """``sn_mesh.bc`` entries (the ``_BoundBoundaryOperator`` shim) forward
        the realized law's role so the mesh-wired BCs are BOUNDARY ops."""
        geom = StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
            bcs=(BC.vacuum, BC.reflective),
        )
        mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=4),))
        quad = Quadrature.gauss_legendre(n_ordinates=4)
        sn = SNMesh(mesh, quad, placeholder_materials(ng=1))
        for face in ("xmin", "xmax"):
            bc = sn.bc[face]
            assert bc.block_role is BlockRole.BOUNDARY, face
            assert isinstance(bc, BoundaryOperator), face
            assert not isinstance(bc, FullOperator), face


# ─────────────────────────────────────────────────────────────────────
# Composer + adjoint role DERIVATION (Wave O step O.2b 4.5)
# ─────────────────────────────────────────────────────────────────────


class TestComposerRoleDerivation:
    """``OperatorSum`` / ``ScaledOperator`` / ``_AdjointOperator`` DERIVE the
    block role from their operands — no hand-stamp (4.5 twin-path retirement).

    The role of a SUM is the union of the blocks its summands touch
    (``BULK ⊔ BOUNDARY = FULL``); scaling and the G-adjoint preserve which
    blocks are touched.
    """

    @staticmethod
    def _ops(sn):
        sig = np.ones((sn.ng, *sn.spatial_shape))
        return (
            StreamingOperator(sn),
            MultiplicationOperator.from_mesh(sig, sn),
            SNBoundaryOperator(sn),
        )

    def test_scaled_preserves_role(self) -> None:
        sn = _slab_mesh()
        L, C, B = self._ops(sn)
        assert (-C).block_role is BlockRole.BULK
        assert isinstance(-B, BoundaryOperator)
        assert (2.0 * L).block_role is BlockRole.FULL

    def test_sum_join_is_union_of_blocks(self) -> None:
        sn = _slab_mesh()
        L, C, B = self._ops(sn)
        # BULK ⊔ BULK = BULK (C + C)
        assert (C + C).block_role is BlockRole.BULK
        # FULL ⊔ BULK = FULL (L + C — the InvertibleOperator, derived not stamped)
        assert (L + C).block_role is BlockRole.FULL
        # BULK ⊔ BOUNDARY = FULL (the discriminating mixed join)
        assert (C - B).block_role is BlockRole.FULL
        # the full within-group loss with the boundary sibling
        assert (L + C - B).block_role is BlockRole.FULL
        assert isinstance(L + C - B, FullOperator)

    def test_adjoint_preserves_role(self) -> None:
        sn = _slab_mesh()
        L, C, B = self._ops(sn)
        # The G-adjoint transposes the 2×2 block matrix → role is preserved.
        assert isinstance(L.H, FullOperator)
        assert L.H.block_role is BlockRole.FULL
        assert isinstance(B.H, BoundaryOperator)
        assert B.H.block_role is BlockRole.BOUNDARY
        assert C.H.block_role is BlockRole.BULK
        # the composite loss adjoint stays FULL
        assert isinstance((L + C - B).H, FullOperator)

    def test_unclassified_operand_propagates_none(self) -> None:
        """A sum with an unclassified (``None``-role) operand is unclassified —
        the derivation does not guess a role it cannot justify."""
        from orpheus.numerics.operator import IdentityOperator

        sn = _slab_mesh()
        L, _, _ = self._ops(sn)
        # IdentityOperator carries no block role (generic primitive).
        assert IdentityOperator().block_role is None
        assert (L + IdentityOperator()).block_role is None
