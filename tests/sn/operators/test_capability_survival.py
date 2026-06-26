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
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    BoundaryOperator,
    BulkOperator,
    FullOperator,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.angular_operator import IncomingSourceOperator
from orpheus.sn.boundary_operator import SNBoundaryOperator
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import InvertibleOperator, StreamingOperator
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
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
    ``solve``; subtracting the boundary operator ``B`` breaks that
    dispatch, so the result is a generic :class:`OperatorSum` whose
    closure law DROPS ``solve`` (no general ``(A + B)⁻¹`` from the
    operands) while KEEPING ``apply``. This pin exercises the
    ``@property``-backed :attr:`SNBoundaryOperator.capabilities` through a
    composition — the existing ``(L + C)`` invertible test does not cover
    the ``−B`` arm (verification plan GAP-2)."""

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

    def test_composite_keeps_apply_drops_solve(self) -> None:
        *_, composite = self._loss_minus_boundary()
        assert CAP_APPLY in composite.capabilities
        assert CAP_SOLVE not in composite.capabilities

    def test_composite_transpose_follows_closure_law(self) -> None:
        # apply_transpose propagates iff every operand has it. This pins the
        # property-backed B's transpose contribution survives the join (so a
        # re-typing that altered B's advertised transpose would red here).
        L, C, B, composite = self._loss_minus_boundary()
        both_have_transpose = all(
            CAP_APPLY_TRANSPOSE in op.capabilities for op in (L + C, B)
        )
        assert (CAP_APPLY_TRANSPOSE in composite.capabilities) == both_have_transpose
