r"""The scheme-honest transpose surface — LD POSITIVE since #310 C2.

Phase 2.5 S0 (#280) minted this file as the LD DEFERRAL pin: the 1-D
reverse walk hand-transposed the DD face chain only, so LD refused loudly
at every surface.  **#310 C2 lifted the deferral** — LD registers its UBLD
``AᵀM⁻¹`` batch VJP (``residual_kernel_batch_transpose``), the reverse
walk's Cartesian arm and the reverse-scan's moment arm carry the
``(…, 2^d)`` cotangent tail, and the same rows that pinned the refusal now
pin the POSITIVE surface (the flip-safe rewrite this file was minted to
force — a flag flipped without its walk reds here).

This file pins:

* the scheme trait ``has_transpose_kernel`` — DERIVED, never declared
  (#310 ruling 2): DD ``True`` (both kernels), LD ``True`` (the batch VJP
  covers its slab-only span), base default ``False``;
* the scheme-honest predicate ``StreamingOperator.is_adjointable`` and its
  a∧b propagation through ``(L+C)`` / ``(L+C−B)`` — now POSITIVE for both
  1-D schemes;
* the eager ``.H`` construction for LD (the former ``MissingAdjoint``
  refusal row, flipped);
* the walk entry's Pattern-4 backstop: a cotangent whose spatial-moment
  tail does not match the scheme's raises a typed ``ValueError`` (the
  silent-broadcast hazard the old guard protected against, now guarded
  SHAPE-wise rather than scheme-wise);
* the ORIENTATION factor (2.5a): a multi-D Cartesian mesh still refuses —
  ``is_adjointable = scheme.has_transpose_kernel ∧
  representation.has_transpose_walk`` and the multi-D reverse walk is the
  #310 C3/C4 arm.  vv Mode-8: function-call asserts only.
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

    def test_dd_derives_true(self) -> None:
        require(
            DiamondDifference.has_transpose_kernel is True,
            "DD must derive has_transpose_kernel=True — it registers BOTH "
            "the Cartesian batch VJP and the curvilinear cell-balance "
            "kernel (#310 C1/C2).",
        )

    def test_ld_derives_true(self) -> None:
        require(
            LinearDiscontinuous.has_transpose_kernel is True,
            "LD must derive has_transpose_kernel=True — it registers the "
            "UBLD AᵀM⁻¹ batch VJP, which alone covers its slab-only "
            "forward span (#310 C2; the covering derivation).",
        )


# ═══════════════════════════════════════════════════════════════════════
# The scheme-honest predicate + a∧b propagation
# ═══════════════════════════════════════════════════════════════════════


class TestPredicateHonesty:
    def test_ld_streaming_is_adjointable(self) -> None:
        sn, _ = _slab(scheme=LinearDiscontinuous())
        require(
            StreamingOperator(sn).is_adjointable,
            "L on a 1-D LD mesh must advertise the adjoint axis — the "
            "batch VJP is registered and CumprodScan reverses (#310 C2).",
        )

    def test_dd_streaming_stays_adjointable(self) -> None:
        sn, _ = _slab(scheme=DiamondDifference())
        require(
            StreamingOperator(sn).is_adjointable,
            "L on a DD mesh must still advertise the adjoint axis — the "
            "honest predicate must not over-refuse (anti-pattern #11).",
        )

    def test_ld_composites_propagate_the_capability(self) -> None:
        sn, sig_t = _slab(scheme=LinearDiscontinuous())
        lc = _lc(sn, sig_t)
        require(
            lc.is_adjointable,
            "(L+C) on 1-D LD must be adjointable (OperatorSum a∧b law).",
        )
        require(
            (lc - SNBoundaryOperator(sn)).is_adjointable,
            "(L+C−B) on 1-D LD must be adjointable (a∧b through the sum).",
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
    def test_ld_H_constructs(self) -> None:
        r"""The former ``MissingAdjoint`` refusal row, FLIPPED (#310 C2)."""
        sn, sig_t = _slab(scheme=LinearDiscontinuous())
        require(
            StreamingOperator(sn).H is not None,
            "L.H on 1-D LD must construct — the deferral is lifted.",
        )
        require(
            _lc(sn, sig_t).H is not None,
            "(L+C).H on 1-D LD must construct — the deferral is lifted.",
        )

    def test_dd_H_constructs(self) -> None:
        sn, sig_t = _slab(scheme=DiamondDifference())
        require(
            _lc(sn, sig_t).H is not None,
            "(L+C).H on DD must construct — the positive control.",
        )


class TestWalkEntryGuard:
    def test_ld_direct_transpose_runs(self) -> None:
        r"""The former typed-deferral row, FLIPPED: a DIRECT Euclidean
        ``apply_transpose`` on a well-formed moment-tailed cotangent runs
        end-to-end and returns the moment-tailed pullback (#310 C2)."""
        sn, sig_t = _slab(scheme=LinearDiscontinuous())
        rng = np.random.default_rng(310)
        phi = _ld_composite(sn)
        phi.interior.values[:] = rng.standard_normal(phi.interior.values.shape)
        out = _lc(sn, sig_t).apply_transpose(phi)
        require(
            out.interior.values.shape == phi.interior.values.shape,
            "LD (L+C)ᵀφ must return the moment-tailed bulk shape.",
        )
        require(
            bool(np.all(np.isfinite(out.interior.values))),
            "LD (L+C)ᵀφ produced non-finite bulk values.",
        )

    def test_ld_tailless_cotangent_refuses_loudly(self) -> None:
        r"""The Pattern-4 backstop that REPLACES the scheme guard: a
        cotangent missing the scheme's spatial-moment tail must raise the
        typed shape refusal, never broadcast silently through the batch
        VJP (the hazard the S0 deferral guard originally protected)."""
        sn, sig_t = _slab(scheme=LinearDiscontinuous())
        tailless = FullField(
            interior=AngularFlux.zeros_on(sn),          # NO moment tail
            boundary=AngularBoundaryFlux.zeros_on(sn),
        )
        with pytest.raises(ValueError, match="spatial-moment tail"):
            _lc(sn, sig_t).apply_transpose(tailless)

    def test_ld_bare_streaming_transpose_runs_too(self) -> None:
        r"""The σ-free ``streaming_action_transpose`` funnels through the
        same walk — the positive sibling of the former refusal row."""
        sn, _ = _slab(scheme=LinearDiscontinuous())
        out = StreamingOperator(sn).apply_transpose(_ld_composite(sn))
        require(
            bool(np.all(np.isfinite(out.interior.values))),
            "LD Lᵀφ produced non-finite bulk values.",
        )

    def test_dd_transpose_positive_control(self) -> None:
        r"""DD's transpose still runs end-to-end and returns finite values
        (anti-pattern #11 — the surface must not over-refuse)."""
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
        sig_t, _psi = het_operands(sn)[:2]
        with pytest.raises(MissingAdjoint):
            _ = StreamingOperator(sn).H
        with pytest.raises(MissingAdjoint):
            _ = _lc(sn, sig_t).H

    def test_cart2d_direct_transpose_still_raises_typed_deferral(self) -> None:
        r"""The backstop behind the predicate: a DIRECT Euclidean
        ``apply_transpose`` (bypassing ``.H``) still hits the
        representation's loud multi-D deferral raise."""
        sn = cart2d_2g_nonsquare()
        sig_t, psi = het_operands(sn)[:2]
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
