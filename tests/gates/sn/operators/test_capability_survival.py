r"""B4 carve safety net — capability + block-role SURVIVAL through the
property-based ``SNBoundaryOperator`` and the bare-mixin
``IncomingSourceOperator`` (GAP-1b, GAP-2 of the verification plan).

The operator-inverse-algebra carve re-types the ``LinearOperator`` family
(``Protocol[V]`` → ``Protocol[D, C]``) and reworks the solve axis. RC 2 of
the verification plan flags the nastiest re-typing target: the
recursive ``SNBoundaryOperator.is_adjointable``. A re-typing that
silently turned that property into a stale plain attribute, or dropped a
leaf's advertised set, would change the capability closure of a composite
``(L + C − B)`` — invisible today because nothing pinned the composite's
surviving set, only ``(L + C)``.

These are foundation (software-invariant) pins: capability sets and the
``None`` block-role default are set-equality / enum-identity facts, not
convergence/flux/eigenvalue claims. References are closed-form (the
advertised sets are fixed by the operators' definitions).

See ``.claude/plans/archive/issue_226_b4_operator_generics_verification.md``
(GAP-1b, GAP-2) and ``.claude/plans/archive/operator_inverse_algebra_carve.md`` §4.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.transport.boundary_realization import stamp_boundary_role
from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.operator import NotInvertible, ZeroOperator
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.numerics.operator import (
    BlockRole,
    BoundaryOperator,
    BulkOperator,
    FullOperator,
    OperatorSum,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.realized_operators import ZeroMap
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.problem import SNProblem
from orpheus.sn.operators.radial_characteristic import RadialCharacteristicOperator
from orpheus.sn.operators.streaming import StreamingCollisionOperator, StreamingOperator
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.transport.mesh.material_xs_field import MaterialXSField
from orpheus.numerics.space import FunctionSpace
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_transfer import (
    IsotropicFission,
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.operators.scattering import ScatteringOperator
from orpheus.transport.operators.transfer import LegendreMomentTransfer
from tests._harness.predicates import (
    INVERTIBLE,
    STRUCTURAL_ABSENT,
    VALUE_RAISE,
    assert_inverse_adjoint_contract,
)
from tests.gates.sn._test_helpers import (
    rc_march,
    material_xs_from_raw,
    placeholder_materials,
)
from orpheus.numerics.basis.spherical_harmonic_basis import SphericalHarmonicBasis
from orpheus.transport.material_field import TransferMaterialField
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


def _sphere_mesh(nx: int = 4, n_ord: int = 4, ng: int = 1) -> SNProblem:
    # The seed-carrying geometry (R12a): the only mesh on which System B —
    # and therefore A_BB — exists.
    geom = StructuredGeometry(
        coord=CoordSystem.SPHERICAL,
        breakpoints=(0.0, 2.0),
        mat_ids=(0,),
        face_roles=(0,),
    )
    mesh = Mesher(geom).partition(CellsByCount.uniform_volume(nx)).mesh
    quad = Quadrature.gauss_legendre(n_ordinates=n_ord)
    return SNProblem(mesh, quad, placeholder_materials(ng=ng), boundary_conditions=conditions_for(mesh, _CONDITIONS))


# ─────────────────────────────────────────────────────────────────────
# GAP-1b — IncomingSourceOperator (bare-mixin) is unclassified by default
# ─────────────────────────────────────────────────────────────────────


class TestTheBlockRoleIsAnInstanceStampNotAClassFact:
    r"""The bare-mixin ``None`` default survives the ``[D, C]`` re-typing —
    pinned on the class that carries BOTH answers in production (RC 3 ∩ RC 1).

    ⭐ **Why this subject, since P3.** This class previously asked the retired
    ``IncomingSourceOperator`` whether it was unclassified, on the claim that
    the rank-0 affine inflow source is ``q.boundary`` and not a linear boundary
    operator ``B``. P3 refuted that claim: prescribed inflow's realization IS a
    linear boundary operator — the ZERO one — and it is stamped like every
    other law (see ``tests/gates/sn/operators/test_operator_block_role.py``, where
    prescribed now sits in ``_LINEAR_LAWS``).

    What survives the retirement, and is what this module actually needs, is the
    orthogonal claim: :attr:`block_role` is an **instance stamp**, so the same
    class reads ``None`` bare and ``BOUNDARY`` once
    :func:`~orpheus.transport.boundary_realization.stamp_boundary_role` has been applied.
    Post-P3 that is no longer hypothetical — :class:`ZeroOperator` genuinely
    appears both ways in production: unstamped as the within-group zero fission
    slot, and stamped as the realization of vacuum / prescribed inflow / a
    re-emission law at ``α = 0``. A re-typing that turned the stamp into a class
    attribute would collapse the two, and this is the row that sees it.
    """

    def test_an_unstamped_instance_is_unclassified(self) -> None:
        op = ZeroOperator()
        assert op.block_role is None
        assert not isinstance(op, BulkOperator)
        assert not isinstance(op, FullOperator)
        assert not isinstance(op, BoundaryOperator)

    def test_stamping_one_instance_leaves_the_class_alone(self) -> None:
        """The negative leg: the stamp must not leak to a sibling instance.

        Without this, ``block_role`` could be a CLASS attribute set by
        ``stamp_boundary_role`` and every row above would still pass — read in
        the wrong order, because the bare instance would be constructed before
        the stamp ever ran.
        """
        stamped = stamp_boundary_role(ZeroOperator())
        assert stamped.block_role is BlockRole.BOUNDARY
        assert isinstance(stamped, BoundaryOperator)
        # …and a FRESH sibling, built after the stamp, is still bare.
        assert ZeroOperator().block_role is None
        assert not isinstance(ZeroOperator(), BoundaryOperator)


# ─────────────────────────────────────────────────────────────────────
# GAP-2a — per-leaf capability set (strict equality, not membership)
# ─────────────────────────────────────────────────────────────────────


class TestRealizedPrescribedInflowCapabilities:
    r"""RC 2: the realized prescribed-inflow surface, per-leaf and strict.

    ⭐ **The surface INVERTED at P3, and the inversion is the claim.** The
    retired ``IncomingSourceOperator`` was apply-only — ``is_adjointable`` was
    ``False``, because an affine map has no transpose to advertise. The zero
    morphism that replaced it is adjointable (``0ᵀ = 0``), so the honest pin is
    the new set, not a preserved one. What survives verbatim is the *shape* of
    the concern this module exists for: a re-typing that dropped a leaf's
    advertised set, or added one it cannot honour, changes the capability
    closure of the composite ``(L + C − B)``.

    ``is_invertible`` stays ``False`` and ``inverse`` stays absent — the zero
    map is the singular map par excellence, and
    :class:`~orpheus.numerics.operator.ZeroOperator` declines to declare an
    ``inverse()`` at all rather than shipping a raising stub.
    """

    def _realized(self):
        """The realized operator of a face declaring a pure inflow (a wall with no return), read off the posing path: a ``ZeroMap``
        whose inflow is the binding's source, not the operator's."""
        mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=0, right=MirrorDeck())).partition(CellsByCount.uniform_volume(4)).mesh
        sn = SNProblem(
            mesh, Quadrature.gauss_legendre(n_ordinates=4), placeholder_materials(ng=1),
            boundary_conditions={0: WallCondition.uniform(2.5)},
        )
        op = sn.bc["xmin"]
        if not isinstance(op, ZeroMap):
            pytest.fail(f"a pure-inflow face must realize to a ZeroMap, got {type(op).__name__}")
        return op

    def test_the_transpose_is_now_advertised(self) -> None:
        op = self._realized()
        assert op.is_adjointable
        assert callable(getattr(op, "apply_transpose", None))

    def test_the_inverse_axis_raises(self) -> None:
        """Meaning change at wave 2 S2b: the retired zero morphism declined to declare an ``inverse()`` at all, while every realized
        boundary operator declares one uniformly (``RealizedBoundary.inverse``, one body) that raises ``NotInvertible`` for a
        projecting image. ``is_invertible`` stays ``False`` and the guard is the value arm of the contract."""
        op = self._realized()
        assert not op.is_invertible
        with pytest.raises(NotInvertible):
            op.inverse()
        assert callable(getattr(op, "apply", None))


# ─────────────────────────────────────────────────────────────────────
# GAP-2b — (L + C − B) composite capability SURVIVAL
# ─────────────────────────────────────────────────────────────────────


class TestLossMinusBoundaryCompositeCapabilities:
    r"""The within-group loss with its boundary sibling, ``L + C − B``.

    The fused :class:`StreamingCollisionOperator` atom (built explicitly
    since S2; ``+`` never fuses, ruling 52) advertises ``solve`` (the DIRECT
    sweep). Subtracting the boundary operator ``B`` gives a plain
    :class:`OperatorSum`, which has no inverse of its own and no ``solve``
    verb (ruling 50); its inverse is the Green of the NAMED split,
    ``GreenOperator(L+C, B)`` (the boundary-Jacobi iteration, typed). The
    pin is WHICH inverse each object has — the sweep for the atom, none for
    the sum, the Green for the split — plus the ``apply`` survival.
    It still exercises :attr:`SNBoundaryOperator.is_adjointable` through
    a composition — the existing ``(L + C)`` invertible test does not
    cover the ``−B`` arm (verification plan GAP-2)."""

    def _loss_minus_boundary(self):
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        L = StreamingOperator.pose(sn)
        C = MultiplicationOperator.from_mesh(sigma_t, sn)
        B = SNBoundaryOperator(sn)
        # S2 (ruling 52): the fused atom is built explicitly, so the
        # composite is (L+C) − B over the sweep-invertible atom, as before.
        return L, C, B, StreamingCollisionOperator(L, C) - B

    def test_boundary_operator_exposes_apply(self) -> None:
        # Precondition for the composite to construct (RC 2 — the eager
        # apply-guard must accept the −B arm at composition time).
        _, _, B, _ = self._loss_minus_boundary()
        assert callable(getattr(B, "apply", None))

    def test_l_plus_c_is_invertible(self) -> None:
        # The contrast case: the (L + C) ATOM is sweep-invertible; since S2
        # (ruling 52) the plain sum ``L + C`` is not — the atom is explicit.
        L, C, _, _ = self._loss_minus_boundary()
        lc = StreamingCollisionOperator(L, C)
        assert lc.is_invertible is True
        assert (L + C).is_invertible is False

    def test_composite_has_no_inverse_its_split_has_a_green(self) -> None:
        # INVERTED at S2 (ruling 50): a sum has no inverse of its own (the
        # leading-term rule is gone) and no solve verb; the Green of the
        # NAMED split M = (L+C), N = B inverts it.
        from orpheus.numerics.green_operator import GreenOperator
        from orpheus.numerics.operator import NotInvertible, OperatorSum

        L, C, B, composite = self._loss_minus_boundary()
        assert composite.is_invertible is False
        with pytest.raises(NotInvertible, match="GreenOperator\\(M, N\\)"):
            composite.inverse()
        assert "solve" not in OperatorSum.__dict__
        lc = next(w for _, w in composite.terms if isinstance(w, StreamingCollisionOperator))
        green = GreenOperator(lc, B)
        assert green.inverse() == composite

    def test_composite_transpose_follows_closure_law(self) -> None:
        # is_adjointable propagates iff every operand is. This pins that
        # B's recursive-predicate contribution survives the join (so a
        # re-typing that altered B's adjointability would red here).
        L, C, B, composite = self._loss_minus_boundary()
        both_adjointable = StreamingCollisionOperator(L, C).is_adjointable and B.is_adjointable
        assert composite.is_adjointable == both_adjointable


# ─────────────────────────────────────────────────────────────────────
# Predicate FAITHFULNESS — the carve keystone (verification spec §2.3)
# ─────────────────────────────────────────────────────────────────────
#
# Phase 2b (rewired at carve P4): the per-operator two-axis CONTRACT for
# the transport + SN advertisers — keystone v2 (spec §36). The numerics
# leaves/composers are pinned in
# ``tests/gates/numerics/test_operator_capability_predicates.py``; THIS file
# closes the transport energy operators + the SN streaming/boundary
# family. Both share the ONE contract assertion
# ``tests/_harness/predicates.assert_inverse_adjoint_contract``.

_SIGS0 = np.array([[0.20, 0.00], [0.05, 0.18]])   # P0 group-transfer (asymmetric)
_SIGS1 = np.array([[0.02, 0.00], [0.01, 0.015]])  # P1 (small anisotropy)
_SIG2 = np.array([[0.00, 0.03], [0.01, 0.00]])    # (n,2n)


def _synthetic_mat_xs(nx: int = 4) -> MaterialXSField:
    """Single-material 2G synthetic XS field (asymmetric SigS, nonzero Sig2)."""
    cells = {0: (np.arange(nx), np.zeros(nx, dtype=int))}
    return material_xs_from_raw(
        sig_s={0: [_SIGS0, _SIGS1]}, sig2={0: _SIG2},
        cells_by_mat=cells, ng=2, nx=nx, ny=1,
    )


class TestPredicateFaithfulness:
    r"""The two-axis inverse/adjoint contract for EVERY transport + SN
    advertiser — keystone v2 (spec §36), the permanent successor of the
    frozenset-coexistence scaffold. Spans the ``(invertible × adjointable)``
    quadrants, including the VALUE-dependent asymmetry leaf
    (``MultiplicationOperator(true-zero-coeff)``: adjointable but NOT
    invertible) that breaks a buggy predicate which merely mirrors the
    other axis (spec §0.6/§8)."""

    def _sn_operators(self, sn):
        spatial = (sn.ng, *sn.spatial_shape)
        sigma_t = np.ones(spatial)
        sigma_singular = sigma_t.copy()
        sigma_singular[0, 0] = 0.0  # a TRUE zero → C is singular (min|f| = 0)
        L = StreamingOperator.pose(sn)
        C = MultiplicationOperator.from_mesh(sigma_t, sn)
        C_singular = MultiplicationOperator.from_mesh(sigma_singular, sn)
        B = SNBoundaryOperator(sn)
        ops = [L, C, C_singular, L + C, L + C - B, B]
        # The realized per-face boundary operators (RealizedBoundary):
        # their is_* MUST delegate to the inner operator, not the base default.
        ops += [sn.bc[face] for face in sn.angular_trace.layout.faces]
        return ops

    def _transport_energy_operators(self):
        mat = _synthetic_mat_xs()
        # F's space is MANDATORY (CS4a K2); a hand-built space keeps the
        # capability survey space-agnostic (its subject is predicates,
        # not posing). The iso pair + S are BOUND to the same synthetic
        # space since the S4-amendment: the adjoint leg takes .H, which
        # refuses an unbound non-multiplier (their constructors have been
        # mandatory since CS4c steps 2/3/5 — this comment said "their
        # Optional-space ctor survives until CS4c's mandatory flip" until
        # 2026-09-07; the survey simply never exercises an unbound state
        # on adjointable rows).
        survey_space = FunctionSpace("capability_survey", (2, 4, 1))
        # S's tier-2 mint needs a real posed composite since CS4c step 3
        # (the faces bind the interior); the survey reads predicates
        # only, so a small slab composite serves the row.
        composite = _slab_mesh(ng=2).full_field_space
        return [
            IsotropicScattering.from_material_xs(mat, space=survey_space),
            IsotropicN2N.from_material_xs(mat, space=survey_space),
            # CS4c step 4: the fission ENERGY binding rides the survey
            # space; the ANGULAR binding (frame-conjugated) needs the
            # posed composite, like S.
            IsotropicFission.from_material_xs(mat, space=survey_space),
            FissionOperator.from_solver_data(mat_xs=mat, space=composite),
            LegendreMomentTransfer.on_basis(
            TransferMaterialField.scattering(mat), SphericalHarmonicBasis(L=1), skip_l0=True,
        ),
            LegendreMomentTransfer.on_basis(
            TransferMaterialField.n2n(mat), SphericalHarmonicBasis(L=1), skip_l0=False,
        ),
            ScatteringOperator.from_solver_data(
                mat_xs=mat,
                scattering_order=0,
                space=composite,
            ),
        ]

    def test_axis_separation_is_value_dependent(self) -> None:
        """The two axes must be INDEPENDENTLY correct. A singular
        ``MultiplicationOperator`` is adjointable but NOT invertible; ``L`` is
        the structural twin (adjointable, not invertible); ``(L+C)`` is both —
        a predicate that merely returned the other axis fails on the singular C."""
        sn = _slab_mesh(ng=2)
        spatial = (sn.ng, *sn.spatial_shape)
        sigma_singular = np.ones(spatial)
        sigma_singular[0, 0] = 0.0
        C_singular = MultiplicationOperator.from_mesh(sigma_singular, sn)
        assert C_singular.is_adjointable is True
        assert C_singular.is_invertible is False
        L = StreamingOperator.pose(sn)
        assert L.is_adjointable is True and L.is_invertible is False
        C_ok = MultiplicationOperator.from_mesh(np.ones(spatial), sn)
        lc = StreamingCollisionOperator(L, C_ok)  # S2: the atom (ruling 52)
        assert lc.is_invertible is True
        assert lc.is_adjointable is True

    def test_boundary_law_wrapper_delegates_predicates(self) -> None:
        """A ``RealizedBoundary`` delegates ``is_*`` to its inner operator,
        NOT the base ``LinearOperator`` default. A non-delegating operator
        would report a vacuum/reflective face as ``is_adjointable=False``
        and silently break the ``B`` aggregator's
        ``all(law.is_adjointable …)`` rule."""
        sn = _slab_mesh(ng=2)
        for face in sn.angular_trace.layout.faces:
            law = sn.bc[face]
            assert law.is_adjointable == law.inner.is_adjointable
            assert law.is_invertible == law.inner.is_invertible

    def test_inverse_adjoint_contract_keystone_v2_sn(self) -> None:
        """KEYSTONE v2, SN/transport slice (carve P4, spec §36).

        The permanent successor of the frozenset scaffold above: pins the
        two-axis contract (returns-vs-raises-vs-absent on the inverse
        axis, eager-``.H`` on the adjoint axis, TypeGuard-bridge
        consistency) over the SN + transport advertisers. Expectations
        VERIFIED against the live fixture (probe 2026-07-02): the vacuum
        realized face operators declare the forwarded ``inverse()`` and their
        guard raises ``NotInvertible`` (VALUE arm); ``L`` and ``B`` are
        the structural arm; the sweep-invertible ``(L+C)`` atom and the Green
        of the split ``M = L+C, N = B`` return, while the plain sum
        ``(L+C) − B`` refuses (since S2, ruling 50).
        """
        sn = _slab_mesh(ng=2)
        spatial = (sn.ng, *sn.spatial_shape)
        L = StreamingOperator.pose(sn)
        C = MultiplicationOperator.from_mesh(np.ones(spatial), sn)
        sigma_singular = np.ones(spatial)
        sigma_singular[0, 0] = 0.0
        C_singular = MultiplicationOperator.from_mesh(sigma_singular, sn)
        B = SNBoundaryOperator(sn)
        lc = StreamingCollisionOperator(L, C)  # S2: the atom, built explicitly
        rows = [
            (L, False, True, STRUCTURAL_ABSENT),
            (C, True, True, INVERTIBLE),
            (C_singular, False, True, VALUE_RAISE),
            (lc, True, True, INVERTIBLE),           # sweep-invertible
            (lc - B, False, True, VALUE_RAISE),     # S2: a sum has no inverse of its own
            (GreenOperator(lc, B), True, True, INVERTIBLE),  # the Green of the split
            (B, False, True, STRUCTURAL_ABSENT),
            (lc.inverse(), True, True, INVERTIBLE),   # SweepOperator (#280 2.5c: adj)
        ]
        # The realized face-law wrappers: forwarded inverse() + guard.
        rows += [
            (sn.bc[face], False, True, VALUE_RAISE)
            for face in sn.angular_trace.layout.faces
        ]
        # Transport energy leaves: all structurally non-invertible,
        # all adjointable (the S†/F† axis, #276).
        rows += [
            (op, False, True, STRUCTURAL_ABSENT)
            for op in self._transport_energy_operators()
        ]
        # A_BB — System B's self-block (coupled-block campaign step 4b):
        # exists ONLY on the seed-carrying sphere, where the direct Carlson
        # march IS its exact inverse and both transposes are realized, so
        # the positive surface is (invertible, adjointable) = (True, True).
        # The canonical durable pin the step-4b follow-up named — gate 5c
        # of TestA_BB_Forward catches a predicate flip-back dynamically;
        # this row is the static table entry beside its siblings.
        sph = _sphere_mesh(ng=2)
        sigma_sph = MultiplicationOperator.from_mesh(
            np.ones((sph.ng, *sph.spatial_shape)), sph
        ).coefficient
        rows.append(
            (rc_march(sph, sigma_sph), True, True, INVERTIBLE)
        )
        for op, inv, adj, contract in rows:
            assert_inverse_adjoint_contract(
                op, invertible=inv, adjointable=adj, inverse_contract=contract
            )

    def test_sweep_vs_green_inverse_keyed_by_type(self) -> None:
        """RE-POSED at S2 (rulings 50, 52): the ``(L+C)`` ATOM (built
        explicitly; ``+`` never fuses) owns the DIRECT sweep, its
        ``.inverse()`` a :class:`SweepOperator`. The plain sum ``(L+C) − B``
        has no inverse of its own; the ITERATIVE :class:`GreenOperator` of
        the NAMED split ``M = L+C, N = B`` inverts it. (Until S2 the operand
        spelling selected the algorithm: ``L + C`` fused by dispatch, and the
        sum's inverse was the Green keyed on its leading term.)
        """
        from orpheus.numerics.operator import NotInvertible

        sn = _slab_mesh(ng=2)
        L = StreamingOperator.pose(sn)
        C = MultiplicationOperator.from_mesh(np.ones((sn.ng, *sn.spatial_shape)), sn)
        lc = StreamingCollisionOperator(L, C)
        assert lc.is_invertible is True and lc.is_adjointable is True
        assert isinstance(lc.inverse(), SweepOperator)  # the atom's own sweep
        B = SNBoundaryOperator(sn)
        composite = lc - B
        assert type(composite) is OperatorSum
        with pytest.raises(NotInvertible, match="GreenOperator\\(M, N\\)"):
            composite.inverse()
        green = GreenOperator(lc, B)
        assert green.inverse() is green.inner      # object-identity involution
        assert green.inverse() == composite        # equal to the sum spelled apart
