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
pinned in ``tests/gates/numerics/test_operator_protocols.py``; here we pin the
SN tagging itself.

For ``C``, ``L`` and ``(L+C)`` we construct instances and assert
``isinstance`` + exclusivity. For ``S`` and ``F`` (which need a
``MaterialXSField``) we assert the class-level tag: because ``block_role``
is a class attribute and the marker metaclass reads the *instance's*
``block_role``, ``isinstance(s_instance, BulkOperator)`` is by
construction equivalent to ``ScatteringOperator.block_role is
BlockRole.BULK`` — the class-level assertion fully determines the
instance membership, so it is not a coverage gap.

The BOUNDARY role tags the realized boundary operators
(``problem.bc`` entries, the :class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`
leaves: vacuum / reflective / white / specular wall / periodic carry
:attr:`BlockRole.BOUNDARY`, and so does a pure inflow, whose linear factor is the
zero map: its source is ``q.boundary``, not part of ``B``). Pinned in
``TestBoundaryLeaves`` below on the entries the posing path wires (wave 2 S2b:
the realizer's raw output and the forwarding shim of the earlier spelling are
the same objects now), each row posed through a real one-face-under-test slab. The O.4a.1-γ tagging lands the role; O.4a.2 wires
``B`` into ``(L_full + C − S − F − B)`` as a sibling of ``L``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.operator import (
    BlockRole,
    BoundaryOperator,
    BulkOperator,
    FullOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.problem import SNProblem
from orpheus.transport.operators.fission import FissionOperator
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.operators.streaming import (
    StreamingCollisionOperator,
    StreamingOperator,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from tests.gates.sn._test_helpers import placeholder_materials
from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


pytestmark = [pytest.mark.foundation]


def _slab_mesh(nx: int = 4, n_ord: int = 4, ng: int = 1) -> SNProblem:
    geom = StructuredGeometry(
        coord=CoordSystem.CARTESIAN,
        breakpoints=(0.0, 2.0),
        mat_ids=(0,),
        face_roles=(0, 0),
    )
    mesh = Mesher(geom).partition(CellsByCount.uniform_volume(nx)).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNProblem(mesh, quad, placeholder_materials(ng=ng), boundary_conditions=conditions_for(mesh, _CONDITIONS))


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
        composite = StreamingCollisionOperator(  # S2: the atom built explicitly
            StreamingOperator.pose(sn), MultiplicationOperator.from_mesh(sigma_t, sn),
        )
        assert composite.block_role is BlockRole.FULL
        assert isinstance(composite, FullOperator)
        assert not isinstance(composite, BulkOperator)


# ─────────────────────────────────────────────────────────────────────
# Boundary leaves (Wave O step O.4a.1-γ)
# ─────────────────────────────────────────────────────────────────────


#: The face under test of every row below: ``xmin`` of a one-slab, whose other face is its partner (the same deck for a
#: periodic pair, a vacuum wall otherwise).
_FACE = "xmin"


def _boundary_under_test(declared):
    """The realized operator of ``xmin`` when it declares ``declared`` (a ``WallCondition``, carried on a tag, or a deck), posed
    through a real slab: the periodic pair is declared on BOTH faces (a translation quotient is symmetric); any other face has a
    vacuum partner."""
    from orpheus.numerics.quadrature import Quadrature as _Q

    conditions: dict[int, WallCondition] = {}
    if isinstance(declared, WallCondition):
        conditions[0] = declared
        left, right = 0, 1
        conditions[1] = WallCondition.vacuum()
    elif isinstance(declared, TranslationDeck):
        left = right = declared
    else:
        left, right = declared, 0
        conditions[0] = WallCondition.vacuum()
    mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=left, right=right)).partition(CellsByCount.uniform_volume(4)).mesh
    sn = SNProblem(mesh, _Q.gauss_legendre(n_ordinates=4), placeholder_materials(ng=1), boundary_conditions=conditions)
    return sn.bc[_FACE]


# Every LINEAR boundary condition + a representative amplitude per kind. Each realizes to a BOUNDARY-block operator (``A_ss`` only).
# The amplitude 0 and 1 hit the zero map and the unit ends; the 0.7 row is the only specular row with a scaled permutation.
_LINEAR_LAWS = {
    "vacuum": WallCondition.vacuum(),
    "reflective_albedo1": MirrorDeck(),
    "white_albedo1": WallCondition.white(1.0),
    "white_albedo05": WallCondition.white(0.5),
    # All three specular amplitudes are kept because they are three DIFFERENT realizations (the zero map, a scaled ordinate
    # permutation, the unit one), and this file is where the BLOCK ROLE claim lives.
    "albedo_specular_0": WallCondition.specular(0.0),
    "albedo_specular_07": WallCondition.specular(0.7),
    "albedo_specular_1": WallCondition.specular(1.0),
    "albedo_white_03": WallCondition.white(0.3),
    "periodic": TranslationDeck(),
    # ⭐ P3 (2026-08-05) — prescribed inflow JOINED this table: the affine law is ``γ₋ψ = L γ₊ψ + q``, this tier realizes ``L``,
    # and for a prescribed inflow ``L = 0``, the SAME object vacuum realizes to. The source travels the boundary-source channel.
    "prescribed_inflow": WallCondition.uniform(2.5),
    # …and the sourceless spelling, whose ``q`` is zero: it must realize to the same LINEAR object.
    "prescribed_no_source": WallCondition.uniform(0.0),
}


class TestBoundaryLeaves:
    """The realized boundary laws advertise :attr:`BlockRole.BOUNDARY`."""

    @pytest.mark.parametrize("law_id", list(_LINEAR_LAWS))
    def test_realized_bc_advertises_boundary_operator(self, law_id: str) -> None:
        op = _boundary_under_test(_LINEAR_LAWS[law_id])
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
        mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=0, right=1)).partition(CellsByCount.uniform_volume(4)).mesh
        sn = SNProblem(
            mesh, Quadrature.gauss_legendre(n_ordinates=4), placeholder_materials(ng=1),
            boundary_conditions={0: WallCondition.vacuum(), 1: WallCondition.uniform(2.5)},
        )
        vacuum, prescribed = sn.bc["xmin"], sn.bc["xmax"]
        assert type(prescribed) is type(vacuum)
        # one trace, so the identity of the bound spaces is available: each operator is bound to its OWN face's half-traces.
        for face, op in (("xmin", vacuum), ("xmax", prescribed)):
            assert op.domain is sn.angular_trace.outflow_space(face)
            assert op.codomain is sn.angular_trace.inflow_space(face)
        # …and the source amplitude is invisible here, which is the point: the LINEAR factor cannot depend on q.
        sourceless = _boundary_under_test(WallCondition.uniform(0.0))
        assert type(sourceless) is type(prescribed)

    def test_mesh_bc_forwards_boundary_role(self) -> None:
        """``problem.bc`` entries carry the BOUNDARY role, so the mesh-wired BCs are BOUNDARY ops."""
        geom = StructuredGeometry(
            coord=CoordSystem.CARTESIAN,
            breakpoints=(0.0, 2.0),
            mat_ids=(0,),
            face_roles=(0, MirrorDeck()),
        )
        mesh = Mesher(geom).partition(CellsByCount.uniform_volume(4)).mesh
        quad = Quadrature.gauss_legendre(n_ordinates=4)
        sn = SNProblem(mesh, quad, placeholder_materials(ng=1), boundary_conditions=conditions_for(mesh, _CONDITIONS))
        for face in ("xmin", "xmax"):
            bc = sn.bc[face]
            assert bc.block_role is BlockRole.BOUNDARY, face
            assert isinstance(bc, BoundaryOperator), face
            assert not isinstance(bc, FullOperator), face


# ─────────────────────────────────────────────────────────────────────
# Composer + adjoint role DERIVATION (Wave O step O.2b 4.5)
# ─────────────────────────────────────────────────────────────────────


class TestComposerRoleDerivation:
    """``OperatorSum`` (incl. the one-term ``α·A``) / ``AdjointOperator`` DERIVE the
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
        # BULK ⊔ BULK = BULK (C + C′): two DISTINCT objects — since S2 ``C + C``
        # collects to the one term 2C and exercises no join.
        C2 = MultiplicationOperator.from_mesh(np.ones((sn.ng, *sn.spatial_shape)), sn)
        assert (C + C2).block_role is BlockRole.BULK
        # FULL ⊔ BULK = FULL (L + C, the plain sum since S2; derived not stamped)
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
