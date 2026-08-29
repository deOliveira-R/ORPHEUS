r"""O.1 (Issue #208 / Wave O) — SN operator-leaf block-role tags.

Pins that every SN operator leaf advertises the correct block-role
classification on the 2×2 bulk/boundary block structure:

* Collision ``C``, scattering ``S``, fission ``F`` → :attr:`BlockRole.BULK`
  (``A_bb`` only; read/write the bulk flux, no boundary action).
* Streaming ``L`` → :attr:`BlockRole.FULL` (couples bulk ↔ boundary —
  reads the inflow trace to seed the sweep, writes the outflow trace).
* ``(L + C)`` :class:`StreamingCollisionOperator` → :attr:`BlockRole.FULL`
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
    ConstantInflowSource,
    IsotropicReturn,
    NoSource,
    PeriodicBoundary,
    PrescribedInflow,
    SpecularReturn,
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
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.method_space import SNMethodSpace
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.operators.fission import FissionOperator
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.operators.streaming import (
    StreamingCollisionOperator,
    StreamingOperator,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from tests.sn._test_helpers import face_method_space, placeholder_materials

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
        L = StreamingOperator.pose(sn)
        assert L.block_role is BlockRole.FULL
        assert isinstance(L, FullOperator)
        assert not isinstance(L, BulkOperator)

    def test_invertible_L_plus_C_is_full(self) -> None:
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        composite = StreamingOperator.pose(sn) + MultiplicationOperator.from_mesh(sigma_t, sn)
        assert isinstance(composite, StreamingCollisionOperator)
        assert composite.block_role is BlockRole.FULL
        assert isinstance(composite, FullOperator)
        assert not isinstance(composite, BulkOperator)


# ─────────────────────────────────────────────────────────────────────
# Boundary leaves (Wave O step O.4a.1-γ)
# ─────────────────────────────────────────────────────────────────────


#: The face every law in this module is installed on. Named once because a
#: law that declares its own orientation (white) must AGREE with it — see
#: :data:`_LINEAR_LAWS`.
_FACE = "xmin"          # ⇒ axis "x", outward normal −ê_x ⇒ outward_sign = −1
_FACE_AXIS = "x"
_FACE_OUTWARD_SIGN = -1


def _boundary_method_space(n_ord: int = 4) -> SNMethodSpace:
    r"""A method space carrying BOTH half-traces, on :data:`_FACE`.

    B3.2 migration: a boundary law is typed :math:`\Gamma_+ \to \Gamma_-`, so
    every branch needs the face's **outflow** ordinates (its domain) as well
    as its inflow ordinates (its codomain). Since B3.4a that includes white
    and prescribed inflow; only albedo and periodic still read ``quadrature``
    alone (B3.4b / B3.4c).
    """
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return face_method_space(quad, face=_FACE)


# Every LINEAR boundary law + a representative amplitude per kind. Each
# realizes to a BOUNDARY-block operator (``A_ss`` only). albedo=0/1 hit the
# Zero/Identity fast paths; α∉{0,1} hits the ScaledOperator path — all three
# must carry the role.
_LINEAR_LAWS = {
    "vacuum": VacuumInflow(),
    "reflective_albedo1": ReflectiveBoundary(axis="x", albedo=1.0),
    "reflective_albedo07": ReflectiveBoundary(axis="x", albedo=0.7),
    # ⚠ The orientation MUST match _FACE. Until B3.4a nothing checked it, and
    # these two rows had shipped `outward_sign=+1` (i.e. "xmax") while being
    # installed on "xmin" — so the realized Lambertian averaged the
    # INSTALLATION face's INFLOW hemisphere and reported nothing. B3.4a's
    # cross-check turned that into a loud BoundaryError; sourcing the
    # orientation from the face constant keeps the two from drifting again.
    "white_albedo1": WhiteBoundary(
        axis=_FACE_AXIS, outward_sign=_FACE_OUTWARD_SIGN, albedo=1.0,
    ),
    "white_albedo05": WhiteBoundary(
        axis=_FACE_AXIS, outward_sign=_FACE_OUTWARD_SIGN, albedo=0.5,
    ),
    # ⚠ B3.4b — these three MIGRATED from the bare spelling
    # ``AlbedoBoundary(albedo=α)``, which SN now REFUSES: with ``G = id`` and
    # ``R = α·I`` nothing names a pairing between the two half-traces, so the
    # law is under-determined on an angular trace. The completed spellings
    # realize through the SAME bodies reflective/white do, which is why the
    # role claim carries over unchanged.
    #
    # All three amplitudes are kept because they remain three DIFFERENT
    # production branches — and this file is where the realized TYPE claim
    # lives, so the α=0 row now pins ``ZeroOperator`` (the narrowed zero map
    # B3.4b opened) rather than the retired full-face one.
    "albedo_specular_0": AlbedoBoundary(0.0, SpecularReturn(axis=_FACE_AXIS)),
    "albedo_specular_1": AlbedoBoundary(1.0, SpecularReturn(axis=_FACE_AXIS)),
    # The diffuse closure carries its own orientation, and it MUST match
    # _FACE for the same reason white's must — see the note above; B3.4b
    # parametrised that cross-check over both carriers.
    "albedo_isotropic_03": AlbedoBoundary(
        0.3, IsotropicReturn(axis=_FACE_AXIS, outward_sign=_FACE_OUTWARD_SIGN),
    ),
    "periodic": PeriodicBoundary(),
    # ⭐ P3 (2026-08-05) — prescribed inflow JOINED this table. It sat outside,
    # in a dedicated ``test_prescribed_inflow_is_not_boundary_operator`` row
    # asserting the opposite, on the claim that the rank-0 affine source is
    # ``q.boundary`` and not a linear boundary operator ``B``. That claim was
    # refuted: the affine law is ``γ₋ψ = L γ₊ψ + q``, this tier realizes ``L``,
    # and for prescribed inflow ``L = 0`` — a perfectly ordinary linear
    # boundary operator, and the SAME object vacuum realizes to. The source
    # travels the boundary-source channel instead.
    #
    # ⛔ The old exclusion was not merely a mis-labelling: because
    # ``SNBoundaryOperator._face_laws`` applies every face's law with no role
    # filter, the unstamped AFFINE operator reached ``B`` anyway
    # (``|B(0)| = q``, measured), so the source was delivered through ``B``
    # *and* through the source channel — a double count. Membership here is
    # what makes the two channels one.
    "prescribed_inflow": PrescribedInflow(
        source=ConstantInflowSource(value=2.5),
    ),
    # …and the sourceless spelling, whose ``q`` is zero: it must realize to the
    # same LINEAR object, since the source is not what this tier carries.
    "prescribed_no_source": PrescribedInflow(source=NoSource()),
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

    def test_prescribed_inflow_realizes_the_same_object_vacuum_does(self) -> None:
        r"""⭐ P3: vacuum and prescribed inflow differ ONLY in :math:`q`.

        REPLACES ``test_prescribed_inflow_is_not_boundary_operator``, whose
        claim (no block role; not a ``BoundaryOperator``) P3 refuted — see the
        note on the ``prescribed_inflow`` rows of :data:`_LINEAR_LAWS`, which
        now carry the role claim through the parameterized gate above.

        This row adds what a parameterized membership cannot say: that the two
        laws realize to the same TYPE with the same two spaces. The affine
        form is :math:`\gamma_-\psi = L\,\gamma_+\psi + q`; both laws have
        :math:`L = 0`, and the only difference between them lives in a term
        this tier does not carry. If a future edit gave prescribed inflow its
        own zero-ish operator again, the membership rows would still pass and
        this one would not.
        """
        space = _boundary_method_space()
        vacuum = SNBoundaryRealizer().realize(VacuumInflow(), space)
        prescribed = SNBoundaryRealizer().realize(
            PrescribedInflow(source=ConstantInflowSource(value=2.5)), space,
        )
        assert type(prescribed) is type(vacuum)
        assert prescribed.domain is vacuum.domain
        assert prescribed.codomain is vacuum.codomain
        # …and the source amplitude is invisible here, which is the point: the
        # LINEAR factor cannot depend on q.
        sourceless = SNBoundaryRealizer().realize(
            PrescribedInflow(source=NoSource()), space,
        )
        assert type(sourceless) is type(prescribed)
        assert sourceless.codomain is prescribed.codomain

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
            StreamingOperator.pose(sn),
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
        # FULL ⊔ BULK = FULL (L + C — the StreamingCollisionOperator, derived not stamped)
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
