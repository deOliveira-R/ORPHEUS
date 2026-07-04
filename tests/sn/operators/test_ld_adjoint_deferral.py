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
  the honest surface must not over-refuse).

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
from tests.sn._test_helpers import placeholder_materials
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
        bulk=AngularFlux.zeros_on(sn, spatial_moments=2),
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
            bool(np.all(np.isfinite(out.bulk.values))),
            "DD (L+C)ᵀφ produced non-finite bulk values.",
        )
