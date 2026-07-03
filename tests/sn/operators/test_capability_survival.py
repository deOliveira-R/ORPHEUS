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

See ``.claude/plans/issue_226_b4_operator_generics_verification.md``
(GAP-1b, GAP-2) and ``.claude/plans/operator_inverse_algebra_carve.md`` §4.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.geometry.boundary import ConstantInflowSource, NoSource
from orpheus.numerics.green_operator import GreenOperator
from orpheus.numerics.operator import (
    BoundaryOperator,
    BulkOperator,
    FullOperator,
    OperatorSum,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.angular import IncomingSourceOperator
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import InvertibleOperator, StreamingOperator
from orpheus.sn.operators.sweep_operator import SweepOperator
from orpheus.transport.mesh.material_xs_field import MaterialXSField
from orpheus.transport.operators.fission import FissionOperator
from orpheus.transport.operators.isotropic_scattering import (
    IsotropicN2N,
    IsotropicScattering,
)
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.operators.scattering import (
    LegendreMomentScattering,
    N2NMomentOperator,
    ScatteringOperator,
)
from tests._harness.predicates import (
    INVERTIBLE,
    STRUCTURAL_ABSENT,
    VALUE_RAISE,
    assert_inverse_adjoint_contract,
)
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


# ─────────────────────────────────────────────────────────────────────
# GAP-1b — IncomingSourceOperator (bare-mixin) is unclassified by default
# ─────────────────────────────────────────────────────────────────────


class TestIncomingSourceOperatorIsUnclassified:
    """The rank-0 affine inflow source carries NO block role — it is the
    boundary *source* ``q.boundary``, not a linear boundary operator ``B``.
    Pins the bare-mixin ``None`` default survives the ``[D, C]`` re-typing
    (RC 3 ∩ RC 1)."""

    def test_default_block_role_is_none(self) -> None:
        op = IncomingSourceOperator(NoSource())
        assert op.block_role is None
        assert not isinstance(op, BulkOperator)
        assert not isinstance(op, FullOperator)
        assert not isinstance(op, BoundaryOperator)


# ─────────────────────────────────────────────────────────────────────
# GAP-2a — per-leaf capability set (strict equality, not membership)
# ─────────────────────────────────────────────────────────────────────


class TestIncomingSourceOperatorCapabilities:
    """RC 2: the apply-only surface survives intact (carve P4 rewire of
    the strict caps-equality pin: the "spurious added capability" concern
    becomes "spurious predicate True", still caught)."""

    def test_surface_is_exactly_apply_only(self) -> None:
        op = IncomingSourceOperator(ConstantInflowSource(value=2.5))
        assert not op.is_invertible and not op.is_adjointable
        assert not hasattr(op, "inverse")
        assert callable(getattr(op, "apply", None))


# ─────────────────────────────────────────────────────────────────────
# GAP-2b — (L + C − B) composite capability SURVIVAL
# ─────────────────────────────────────────────────────────────────────


class TestLossMinusBoundaryCompositeCapabilities:
    r"""The within-group loss with its boundary sibling, ``L + C − B``.

    ``L + C`` dispatches to :class:`InvertibleOperator` and advertises
    ``solve`` (the DIRECT sweep); subtracting the boundary operator ``B``
    breaks that dispatch, so the result is a generic
    :class:`OperatorSum` — which, since #226 taxonomy step 4, is
    GREEN-invertible: the leading term (the fused ``L+C``) is invertible,
    so :meth:`inverse` derives the preconditioned-splitting
    :class:`~orpheus.numerics.green_operator.GreenOperator` (the
    boundary-Jacobi iteration, typed); as a generic sum it carries NO
    ``solve`` verb (carve P4 — solving is the inverse OBJECT's apply).
    The pin is WHICH inverse the spelling selects — sweep for the fused
    composite, Green for the generic sum — plus the ``apply`` survival.
    It still exercises :attr:`SNBoundaryOperator.is_adjointable` through
    a composition — the existing ``(L + C)`` invertible test does not
    cover the ``−B`` arm (verification plan GAP-2)."""

    def _loss_minus_boundary(self):
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        L = StreamingOperator(sn)
        C = MultiplicationOperator.from_mesh(sigma_t, sn)
        B = SNBoundaryOperator(sn)
        return L, C, B, L + C - B

    def test_boundary_operator_exposes_apply(self) -> None:
        # Precondition for the composite to construct (RC 2 — the eager
        # apply-guard must accept the −B arm at composition time).
        _, _, B, _ = self._loss_minus_boundary()
        assert callable(getattr(B, "apply", None))

    def test_l_plus_c_is_invertible(self) -> None:
        # The contrast case: (L + C) is sweep-invertible (InvertibleOperator).
        L, C, _, _ = self._loss_minus_boundary()
        lc = L + C
        assert isinstance(lc, InvertibleOperator)
        assert lc.is_invertible is True

    def test_composite_is_green_invertible_without_solve(self) -> None:
        # MIGRATED (#226 step 4; rewired carve P4): the generic sum HAS an
        # inverse — the Green splitting keyed on its invertible leading
        # term — and, as a generic sum, NO solve verb of its own (Design
        # B: the sweep subclass keeps its override; solving with a plain
        # OperatorSum is the inverse OBJECT's apply).
        from orpheus.numerics.green_operator import GreenOperator
        from orpheus.numerics.operator import OperatorSum

        *_, composite = self._loss_minus_boundary()
        assert composite.is_invertible is True
        assert type(composite.inverse()) is GreenOperator
        assert "solve" not in OperatorSum.__dict__

    def test_composite_transpose_follows_closure_law(self) -> None:
        # is_adjointable propagates iff every operand is. This pins that
        # B's recursive-predicate contribution survives the join (so a
        # re-typing that altered B's adjointability would red here).
        L, C, B, composite = self._loss_minus_boundary()
        both_adjointable = (L + C).is_adjointable and B.is_adjointable
        assert composite.is_adjointable == both_adjointable


# ─────────────────────────────────────────────────────────────────────
# Predicate FAITHFULNESS — the carve keystone (verification spec §2.3)
# ─────────────────────────────────────────────────────────────────────
#
# Phase 2b (rewired at carve P4): the per-operator two-axis CONTRACT for
# the transport + SN advertisers — keystone v2 (spec §36). The numerics
# leaves/composers are pinned in
# ``tests/numerics/test_operator_capability_predicates.py``; THIS file
# closes the transport energy operators + the SN streaming/boundary
# family. Both share the ONE contract assertion
# ``tests/_harness/predicates.assert_inverse_adjoint_contract``.

_SIGS0 = np.array([[0.20, 0.00], [0.05, 0.18]])   # P0 group-transfer (asymmetric)
_SIGS1 = np.array([[0.02, 0.00], [0.01, 0.015]])  # P1 (small anisotropy)
_SIG2 = np.array([[0.00, 0.03], [0.01, 0.00]])    # (n,2n)


def _synthetic_mat_xs(nx: int = 4) -> MaterialXSField:
    """Single-material 2G synthetic XS field (asymmetric SigS, nonzero Sig2)."""
    cells = {0: (np.arange(nx), np.zeros(nx, dtype=int))}
    return MaterialXSField._synthetic_for_tests(
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
        L = StreamingOperator(sn)
        C = MultiplicationOperator.from_mesh(sigma_t, sn)
        C_singular = MultiplicationOperator.from_mesh(sigma_singular, sn)
        B = SNBoundaryOperator(sn)
        ops = [L, C, C_singular, L + C, L + C - B, B]
        # The realized per-face boundary-law wrappers (_BoundBoundaryOperator):
        # their is_* MUST delegate to the inner law, not the base default.
        ops += [sn.bc[face] for face in sn.angular_trace.layout.faces]
        return ops

    def _transport_energy_operators(self):
        mat = _synthetic_mat_xs()
        return [
            IsotropicScattering(mat),
            IsotropicN2N(mat),
            FissionOperator(mat_xs=mat),
            LegendreMomentScattering(mat_xs=mat, L=1, skip_l0=True),
            N2NMomentOperator(mat_xs=mat, L=1),
            ScatteringOperator(
                mat_xs=mat,
                quadrature=Quadrature.gauss_legendre(n_ordinates=4),
                scattering_order=0,
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
        L = StreamingOperator(sn)
        assert L.is_adjointable is True and L.is_invertible is False
        C_ok = MultiplicationOperator.from_mesh(np.ones(spatial), sn)
        assert (L + C_ok).is_invertible is True
        assert (L + C_ok).is_adjointable is True

    def test_boundary_law_wrapper_delegates_predicates(self) -> None:
        """The ``_BoundBoundaryOperator`` wrapper delegates ``is_*`` to its inner
        realized law, NOT the base ``LinearOperator`` default. A
        non-delegating wrapper would report a vacuum/reflective face law
        as ``is_adjointable=False`` and silently break the ``B``
        aggregator's ``all(law.is_adjointable …)`` rule."""
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
        face-law shims declare the forwarded ``inverse()`` and their
        guard raises ``NotInvertible`` (VALUE arm); ``L`` and ``B`` are
        the structural arm; the sweep-invertible ``(L+C)`` and the plain
        Green-invertible ``(L+C-B)`` both return.
        """
        sn = _slab_mesh(ng=2)
        spatial = (sn.ng, *sn.spatial_shape)
        L = StreamingOperator(sn)
        C = MultiplicationOperator.from_mesh(np.ones(spatial), sn)
        sigma_singular = np.ones(spatial)
        sigma_singular[0, 0] = 0.0
        C_singular = MultiplicationOperator.from_mesh(sigma_singular, sn)
        B = SNBoundaryOperator(sn)
        rows = [
            (L, False, True, STRUCTURAL_ABSENT),
            (C, True, True, INVERTIBLE),
            (C_singular, False, True, VALUE_RAISE),
            (L + C, True, True, INVERTIBLE),        # sweep-invertible
            (L + C - B, True, True, INVERTIBLE),    # plain sum → Green
            (B, False, True, STRUCTURAL_ABSENT),
            ((L + C).inverse(), True, False, INVERTIBLE),  # SweepOperator
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
        for op, inv, adj, contract in rows:
            assert_inverse_adjoint_contract(
                op, invertible=inv, adjointable=adj, inverse_contract=contract
            )

    def test_sweep_vs_green_inverse_keyed_by_type(self) -> None:
        """``(L+C)`` is the ONE sweep-invertible OperatorSum — its
        ``.inverse()`` override (→ direct :class:`SweepOperator`) shadows
        the generic sum's by MRO.  Subtracting ``B`` breaks the fused
        dispatch: the composite is a PLAIN sum whose inverse is the
        ITERATIVE :class:`GreenOperator` over the sweep-preconditioned
        splitting (#226 step 4's ordering ruling — the operand spelling
        selects the algorithm, and the fused type wins where it exists).

        MIGRATED from ``test_invertible_operator_is_the_sole_invertible_
        sum`` (pre-step-4: ``composite.is_invertible is False`` — "no
        general (A+B)⁻¹").  The sum now HAS a general inverse; what
        ``(L+C)`` remains sole owner of is the DIRECT sweep.
        """
        sn = _slab_mesh(ng=2)
        L = StreamingOperator(sn)
        C = MultiplicationOperator.from_mesh(np.ones((sn.ng, *sn.spatial_shape)), sn)
        lc = L + C
        assert isinstance(lc, InvertibleOperator)
        assert lc.is_invertible is True and lc.is_adjointable is True
        assert isinstance(lc.inverse(), SweepOperator)  # the MRO shadow
        composite = lc - SNBoundaryOperator(sn)
        assert type(composite) is OperatorSum  # the fused dispatch broke
        assert composite.is_invertible is True  # leading term (L+C) invertible
        green = composite.inverse()
        assert isinstance(green, GreenOperator)  # the generic arm → Green
        assert green.inverse() is composite  # object-identity involution
