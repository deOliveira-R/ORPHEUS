r"""B4 carve safety net — capability + block-role SURVIVAL through the
property-based ``SNBoundaryOperator`` and the bare-mixin
``IncomingSourceOperator`` (GAP-1b, GAP-2 of the verification plan).

The operator-inverse-algebra carve re-types the ``LinearOperator`` family
(``Protocol[V]`` → ``Protocol[D, C]``) and reworks the solve axis. RC 2 of
the verification plan flags the nastiest re-typing target: the
``@property``-backed ``SNBoundaryOperator.capabilities``. A re-typing that
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
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
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
from tests._harness.predicates import assert_capability_faithful
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
    """RC 2: the ``ClassVar[frozenset]`` capability set survives intact.

    ``test_bc_universal_invariants.py`` pins this by membership; the strict
    equality below additionally catches a re-typing that ADDS a spurious
    capability tag (the membership test only catches a DROP of apply or a
    spurious solve/transpose)."""

    def test_capabilities_are_exactly_apply(self) -> None:
        op = IncomingSourceOperator(ConstantInflowSource(value=2.5))
        assert op.capabilities == frozenset({CAP_APPLY})


# ─────────────────────────────────────────────────────────────────────
# GAP-2b — (L + C − B) composite capability SURVIVAL
# ─────────────────────────────────────────────────────────────────────


class TestLossMinusBoundaryCompositeCapabilities:
    r"""The within-group loss with its boundary sibling, ``L + C − B``.

    ``L + C`` dispatches to :class:`InvertibleOperator` and advertises
    ``solve`` (the DIRECT sweep); subtracting the boundary operator ``B``
    breaks that dispatch, so the result is a generic
    :class:`OperatorSum` — which, since #226 taxonomy step 4, carries its
    OWN ``solve``: the leading term (the fused ``L+C``) is invertible, so
    the closure derives the preconditioned-splitting
    :class:`~orpheus.numerics.green_operator.GreenOperator` inverse (the
    boundary-Jacobi iteration, typed).  The pin is now WHICH inverse the
    spelling selects — sweep for the fused composite, Green for the
    generic sum — plus the ``apply`` survival.  It still exercises the
    ``@property``-backed :attr:`SNBoundaryOperator.capabilities` through
    a composition — the existing ``(L + C)`` invertible test does not
    cover the ``−B`` arm (verification plan GAP-2)."""

    def _loss_minus_boundary(self):
        sn = _slab_mesh()
        sigma_t = np.ones((sn.ng, *sn.spatial_shape))
        L = StreamingOperator(sn)
        C = MultiplicationOperator.from_mesh(sigma_t, sn)
        B = SNBoundaryOperator(sn)
        return L, C, B, L + C - B

    def test_boundary_operator_advertises_apply(self) -> None:
        # Precondition for the composite to construct (RC 2 — the property
        # set must carry apply for the OperatorSum to accept the −B arm).
        _, _, B, _ = self._loss_minus_boundary()
        assert CAP_APPLY in B.capabilities

    def test_l_plus_c_is_invertible_with_solve(self) -> None:
        # The contrast case: (L + C) HAS solve via InvertibleOperator.
        L, C, _, _ = self._loss_minus_boundary()
        lc = L + C
        assert isinstance(lc, InvertibleOperator)
        assert CAP_SOLVE in lc.capabilities

    def test_composite_keeps_apply_and_carries_green_solve(self) -> None:
        # MIGRATED (#226 step 4): pre-step-4 this pinned "CAP_SOLVE not in
        # caps — no general (A+B)⁻¹".  The generic sum now HAS an inverse
        # — the Green splitting keyed on its invertible leading term — so
        # CAP_SOLVE survives the −B arm (faithfulness: is_invertible ≡
        # CAP_SOLVE).
        *_, composite = self._loss_minus_boundary()
        assert CAP_APPLY in composite.capabilities
        assert CAP_SOLVE in composite.capabilities
        assert composite.is_invertible is True

    def test_composite_transpose_follows_closure_law(self) -> None:
        # apply_transpose propagates iff every operand has it. This pins the
        # property-backed B's transpose contribution survives the join (so a
        # re-typing that altered B's advertised transpose would red here).
        L, C, B, composite = self._loss_minus_boundary()
        both_have_transpose = all(
            CAP_APPLY_TRANSPOSE in op.capabilities for op in (L + C, B)
        )
        assert (CAP_APPLY_TRANSPOSE in composite.capabilities) == both_have_transpose


# ─────────────────────────────────────────────────────────────────────
# Predicate FAITHFULNESS — the carve keystone (verification spec §2.3)
# ─────────────────────────────────────────────────────────────────────
#
# Phase 2b extends the capability-SURVIVAL net above with the per-operator
# predicate-FAITHFULNESS invariant ``is_X ≡ CAP_X ∈ capabilities`` for the
# transport + SN advertisers. This is the equivalence that licenses deleting
# the frozenset in Phase 4: the numerics leaves/composers are pinned in
# ``tests/numerics/test_operator_capability_predicates.py`` and the frame faces
# in ``tests/numerics/test_frame.py``; THIS file closes the transport energy
# operators + the SN streaming/boundary family. All three share the ONE
# keystone assertion ``tests/_harness/predicates.assert_capability_faithful``.

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
    r"""``is_invertible ≡ CAP_SOLVE∈caps`` AND ``is_adjointable ≡
    CAP_APPLY_TRANSPOSE∈caps`` for EVERY transport + SN advertiser — the
    keystone (verification spec §2.3) that licenses Phase 4's frozenset
    deletion. Spans the ``(invertible × adjointable)`` quadrants, including the
    VALUE-dependent asymmetry leaf (``MultiplicationOperator(true-zero-coeff)``:
    adjointable but NOT invertible) that breaks a buggy predicate which merely
    mirrors the other axis (spec §0.6/§8)."""

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
        ops += [sn.bc[face] for face in sn.trace.layout.faces]
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

    def test_sn_operators_are_faithful(self) -> None:
        for op in self._sn_operators(_slab_mesh(ng=2)):
            assert_capability_faithful(op)

    def test_transport_energy_operators_are_faithful(self) -> None:
        for op in self._transport_energy_operators():
            assert_capability_faithful(op)

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
        realized law (mirroring its ``capabilities`` delegation), NOT the base
        ``LinearOperator`` default. A non-delegating wrapper would report a
        vacuum/reflective face law as ``is_adjointable=False`` and silently break
        the ``B`` aggregator's ``all(law.is_adjointable …)`` rule."""
        sn = _slab_mesh(ng=2)
        for face in sn.trace.layout.faces:
            law = sn.bc[face]
            assert law.is_adjointable == (CAP_APPLY_TRANSPOSE in law.capabilities)
            assert law.is_invertible == (CAP_SOLVE in law.capabilities)

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
