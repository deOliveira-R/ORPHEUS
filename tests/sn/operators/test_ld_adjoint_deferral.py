r"""LD adjoint faces are a TYPED deferral — the scheme-honest transpose surface.

Phase 2.5 S0 (#280, FLAG-2 of the 2.5-P0 recon): the 1-D reverse walk
(`_OneDimScanWalk.loss_action_transpose`) hand-transposes the
Diamond-Difference face-flux chain ONLY — its buffers carry no
spatial-moment tail — yet nothing guarded a Linear-Discontinuous slab out
of it: ``CumprodScan.supports`` admits LD slab (``is_affine_scannable``),
``StreamingOperator.is_adjointable`` was unconditionally ``True``, and an
LD ``.H``/``apply_transpose`` either broadcast-crashed or silently
mis-computed on shape coincidence — violating the Protocol's "never a
silent wrong answer" deferral contract.

This file pins the honest surface landed in S0:

* the scheme trait ``has_transpose_kernel`` (DD ``True``, LD ``False``,
  base default ``False`` — a scheme is not assumed transposable until a
  reverse-mode realization of its cell relation exists);
* the scheme-honest predicate ``StreamingOperator.is_adjointable`` and its
  a∧b propagation through ``(L+C)`` / ``(L+C−B)``;
* the eager ``.H`` refusal (``MissingAdjoint`` at CONSTRUCTION — Pattern 4,
  illegal state unrepresentable);
* the reverse walk's entry guard: a DIRECT Euclidean ``apply_transpose``
  call (bypassing ``.H``) raises the typed deferral, never a silent
  broadcast;
* the DD positive control (vv anti-pattern #11: paired positive+negative —
  the honest surface must not over-refuse);
* the ORIENTATION factor (2.5a — the S0-deferred multi-D question):
  ``is_adjointable = scheme.has_transpose_kernel ∧
  representation.has_transpose_walk``, so a multi-D Cartesian mesh (whose
  walks carry no reverse traversal — the #280 deferral) refuses the eager
  ``.H`` at construction exactly like LD, with the representations' loud
  raises as the direct-call backstop.

The deferral lifts with the #280 kernel-pair registration (DD registers
forward + transpose kernels; LD forward-only → the deferral becomes
structural instead of guarded). vv Mode-8: function-call asserts only.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.operator import MissingAdjoint
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.loss_representation import (
    CumprodScan,
    MovingFrontierWindow,
    ScanMarch,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.boundary import SNBoundaryOperator
from orpheus.sn.operators.streaming import StreamingOperator
from orpheus.transport.fields.angular_boundary_flux import AngularBoundaryFlux
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.full_field import FullField
from orpheus.transport.operators.multiplication_operator import MultiplicationOperator
from orpheus.transport.spatial.diamond import DiamondDifference
from orpheus.transport.spatial.linear_discontinuous import LinearDiscontinuous
from orpheus.transport.spatial.scheme import DiscretizationSchemeBase
from tests.sn._test_helpers import (
    cart2d_2g_nonsquare,
    het_operands,
    placeholder_materials,
)
from tests.sn.operators.test_g_adjoint_reciprocity import _random_composite

pytestmark = pytest.mark.foundation


def require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def _slab(scheme=None, ng: int = 2):
    """Reflective 2G slab; per-group-varying σ_t (the reciprocity-file recipe)."""
    quad = Quadrature.gauss_legendre(4)
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5),
        mat_ids=np.zeros(4, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    sn = SNMesh(mesh, quad, placeholder_materials(ng=ng), scheme=scheme)
    sig_t = np.stack(
        [np.full(sn.spatial_shape, 0.5 * (1.0 + 0.5 * g)) for g in range(ng)], axis=0
    )
    return sn, sig_t


def _lc(sn, sig_t):
    return StreamingOperator(sn) + MultiplicationOperator.from_mesh(sig_t, sn)


def _ld_composite(sn) -> FullField:
    """A well-formed LD cotangent: moment-tailed bulk + scalar face trace."""
    return FullField(
        interior=AngularFlux.zeros_on(sn, spatial_moments=2),
        boundary=AngularBoundaryFlux.zeros_on(sn),
    )


# ═══════════════════════════════════════════════════════════════════════
# The scheme trait — declarations pinned
# ═══════════════════════════════════════════════════════════════════════


class TestSchemeTraitDeclarations:
    def test_base_default_is_false(self) -> None:
        require(
            DiscretizationSchemeBase.has_transpose_kernel is False,
            "has_transpose_kernel must default False on the base — a scheme "
            "is NOT assumed transposable until a reverse-mode realization "
            "of its cell relation exists (opt-in, like is_affine_scannable).",
        )

    def test_dd_declares_true(self) -> None:
        require(
            DiamondDifference.has_transpose_kernel is True,
            "DD must declare has_transpose_kernel=True — the 1-D reverse "
            "walk hand-transposes the diamond face-flux chain (O.2b).",
        )

    def test_ld_declares_false(self) -> None:
        require(
            LinearDiscontinuous.has_transpose_kernel is False,
            "LD must declare has_transpose_kernel=False — the UBLD "
            "Schur-residual VJP is the #280 kernel-pair deferral.",
        )


# ═══════════════════════════════════════════════════════════════════════
# The scheme-honest predicate + a∧b propagation
# ═══════════════════════════════════════════════════════════════════════


class TestPredicateHonesty:
    def test_ld_streaming_is_not_adjointable(self) -> None:
        sn, _ = _slab(scheme=LinearDiscontinuous())
        require(
            not StreamingOperator(sn).is_adjointable,
            "L on an LD mesh must NOT advertise the adjoint axis — the "
            "reverse walk carries no LD transpose (names need invariants).",
        )

    def test_dd_streaming_stays_adjointable(self) -> None:
        sn, _ = _slab(scheme=DiamondDifference())
        require(
            StreamingOperator(sn).is_adjointable,
            "L on a DD mesh must still advertise the adjoint axis — the "
            "honest predicate must not over-refuse (anti-pattern #11).",
        )

    def test_ld_composites_propagate_the_refusal(self) -> None:
        sn, sig_t = _slab(scheme=LinearDiscontinuous())
        lc = _lc(sn, sig_t)
        require(
            not lc.is_adjointable,
            "(L+C) on LD must not be adjointable (OperatorSum a∧b law).",
        )
        require(
            not (lc - SNBoundaryOperator(sn)).is_adjointable,
            "(L+C−B) on LD must not be adjointable (a∧b through the sum).",
        )

    def test_dd_composites_stay_adjointable(self) -> None:
        sn, sig_t = _slab(scheme=DiamondDifference())
        lc = _lc(sn, sig_t)
        require(lc.is_adjointable, "(L+C) on DD must stay adjointable.")
        require(
            (lc - SNBoundaryOperator(sn)).is_adjointable,
            "(L+C−B) on DD must stay adjointable.",
        )


# ═══════════════════════════════════════════════════════════════════════
# Eager .H refusal (front door) + the walk entry guard (backstop)
# ═══════════════════════════════════════════════════════════════════════


class TestEagerRefusal:
    def test_ld_H_raises_missing_adjoint_at_construction(self) -> None:
        sn, sig_t = _slab(scheme=LinearDiscontinuous())
        with pytest.raises(MissingAdjoint):
            _ = StreamingOperator(sn).H
        with pytest.raises(MissingAdjoint):
            _ = _lc(sn, sig_t).H

    def test_dd_H_constructs(self) -> None:
        sn, sig_t = _slab(scheme=DiamondDifference())
        require(
            _lc(sn, sig_t).H is not None,
            "(L+C).H on DD must construct — the positive control.",
        )


class TestWalkEntryGuard:
    def test_ld_direct_transpose_raises_typed_deferral(self) -> None:
        r"""A DIRECT Euclidean ``apply_transpose`` (no ``.H``) must refuse
        LOUDLY at the walk entry — the backstop behind the predicate."""
        sn, sig_t = _slab(scheme=LinearDiscontinuous())
        lc = _lc(sn, sig_t)
        with pytest.raises(NotImplementedError, match="no transpose kernel"):
            lc.apply_transpose(_ld_composite(sn))

    def test_ld_bare_streaming_transpose_raises_too(self) -> None:
        r"""The σ-free ``streaming_action_transpose`` funnels through the
        same walk entry — one guard point covers both routes."""
        sn, _ = _slab(scheme=LinearDiscontinuous())
        with pytest.raises(NotImplementedError, match="no transpose kernel"):
            StreamingOperator(sn).apply_transpose(_ld_composite(sn))

    def test_dd_transpose_positive_control(self) -> None:
        r"""DD's transpose still runs end-to-end and returns finite values
        (anti-pattern #11 — the guard must not over-refuse)."""
        sn, sig_t = _slab(scheme=DiamondDifference())
        rng = np.random.default_rng(2026)
        out = _lc(sn, sig_t).apply_transpose(_random_composite(sn, rng))
        require(
            bool(np.all(np.isfinite(out.interior.values))),
            "DD (L+C)ᵀφ produced non-finite bulk values.",
        )


# ═══════════════════════════════════════════════════════════════════════
# The ORIENTATION factor — multi-D honesty (2.5a; the S0-deferred question)
# ═══════════════════════════════════════════════════════════════════════


class TestMultiDOrientationHonesty:
    r"""``is_adjointable`` factorizes along the #280 axes: KERNEL
    (``scheme.has_transpose_kernel``) ∧ ORIENTATION
    (``representation.has_transpose_walk``).  A DD multi-D Cartesian mesh
    passes the kernel factor but its walks carry no reverse traversal —
    pre-2.5a the predicate LIED (``True`` while ``apply_transpose``
    raised); now the eager ``.H`` refuses at construction and the
    representations' raises are the direct-call backstop."""

    def test_rep_trait_declarations(self) -> None:
        sn1, _ = _slab(scheme=DiamondDifference())
        require(
            CumprodScan(sn1).has_transpose_walk is True,
            "CumprodScan (1-D only) must declare has_transpose_walk=True — "
            "the shared 1-D loop walk reverses (2.5a).",
        )
        require(
            ScanMarch(sn1).has_transpose_walk is True,
            "ScanMarch on a 1-D mesh must declare has_transpose_walk=True "
            "— its 1-D branch delegates to the same reverse loop walk.",
        )
        sn2 = cart2d_2g_nonsquare()
        require(
            ScanMarch(sn2).has_transpose_walk is False,
            "ScanMarch on a 2-D mesh must declare has_transpose_walk=False "
            "— the multi-D reverse walk is the #280 deferral.",
        )
        require(
            MovingFrontierWindow(sn2).has_transpose_walk is False,
            "The DAG-wavefront family must inherit has_transpose_walk="
            "False (base opt-in default) — its adjoint is deferred.",
        )

    def test_cart2d_dd_streaming_is_not_adjointable(self) -> None:
        sn = cart2d_2g_nonsquare()
        require(
            not StreamingOperator(sn).is_adjointable,
            "L on a 2-D Cartesian DD mesh must NOT advertise the adjoint "
            "axis — the kernel factor passes but the orientation factor "
            "fails (no multi-D reverse walk); a True here is the predicate "
            "lie the S0 STATUS deferred to 2.5a.",
        )

    def test_cart2d_H_raises_missing_adjoint_at_construction(self) -> None:
        sn = cart2d_2g_nonsquare()
        sig_t, _psi = het_operands(sn)
        with pytest.raises(MissingAdjoint):
            _ = StreamingOperator(sn).H
        with pytest.raises(MissingAdjoint):
            _ = _lc(sn, sig_t).H

    def test_cart2d_direct_transpose_still_raises_typed_deferral(self) -> None:
        r"""The backstop behind the predicate: a DIRECT Euclidean
        ``apply_transpose`` (bypassing ``.H``) still hits the
        representation's loud multi-D deferral raise."""
        sn = cart2d_2g_nonsquare()
        sig_t, psi = het_operands(sn)
        with pytest.raises(NotImplementedError, match="multi-D Cartesian adjoint"):
            _lc(sn, sig_t).apply_transpose(psi)

    def test_1d_dd_positive_control_unmoved(self) -> None:
        r"""The orientation factor must not over-refuse (anti-pattern #11):
        1-D DD keeps the full adjoint surface."""
        sn, sig_t = _slab(scheme=DiamondDifference())
        require(
            StreamingOperator(sn).is_adjointable,
            "1-D DD L must stay adjointable under the two-factor predicate.",
        )
        require(
            _lc(sn, sig_t).H is not None,
            "1-D DD (L+C).H must still construct.",
        )
