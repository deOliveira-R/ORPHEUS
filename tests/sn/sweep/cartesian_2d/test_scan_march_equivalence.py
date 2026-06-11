r"""ScanMarch ≡ oracle — the scan-march schedule (issue #222) principled-equivalence.

The 2-D scan-march sweep (:func:`orpheus.sn.sweep._sweep_2d_scanmarch`, the
:class:`~orpheus.sn.loss_representation.ScanMarch` strategy) is a DIFFERENT
topological linearization of the SAME lower-triangular ``(L+C)`` solve than the
anti-diagonal wavefront: *scan along x, march over y*.  So it is
principled-equivalent (NOT bit-identical) to the
:class:`~orpheus.sn.loss_representation.FullFieldWavefront` ORACLE — they agree to
FP-association.

``vv-principles`` "Bit-identity vs principled-equivalence": the row-march and
the anti-diagonal process the same cell dependencies in a *different order*, so
the converged values differ only by FP non-associativity.  The oracle is the
unconditionally-stable reference (the SAME ``DiamondDifference`` cell kernel via
a different schedule), transitively pinned to the analytical ``k_inf = 1.875``
and ``φ = Q/Σ_t`` grounds.  This is the issue #222 / ``scan_march_verification.md``
G2 gate (the ``FullFieldWavefront`` oracle pins the scan-march at nulp).

Config (Cardinal Rule ≥2G + heterogeneous): a random per-ordinate source
(genuinely anisotropic), heterogeneous ``Σ_t``, random non-zero boundary
inflow, and NON-SQUARE meshes (the x↔y swap moat — Mode-2).  ``foundation`` — a
software invariant, no theory ``:label:``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, Mesh2D
from orpheus.numerics.projection import MomentProjection
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import StreamingOperator
from orpheus.sn.loss_representation import (
    FullFieldWavefront,
    MovingFrontierWindow,
    ScanMarch,
)
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.boundary_flux import BoundaryFlux
from orpheus.transport.timed_full_field import TimedFullField

pytestmark = pytest.mark.foundation

# (nx, ny, level, ng, bc)
CASES = [
    (8, 8, 4, 2, "vacuum"),
    (8, 8, 4, 2, "reflective"),
    (12, 7, 6, 2, "reflective"),    # non-square — the x↔y swap moat (Mode-2)
    (5, 9, 4, 4, "reflective"),     # non-square, 4g
]

# Principled-equivalence tolerance (NOT bit-identity — ``vv-principles``): the
# row-march and the anti-diagonal differ by FP-association BY CONSTRUCTION.  The
# error is bounded ABSOLUTELY at ``~(nx + ny)·eps ≈ 1e-14`` (the x-scan depth +
# the y-march that threads ψ_y row-to-row), so ``assert_allclose(rtol, atol)``
# is the robust metric here — NOT ``assert_array_almost_equal_nulp``, which
# amplifies near-zero boundary elements (a ~1e-15 abs diff at a ~1e-4 shed value
# reads as ~16000 ULP-at-that-scale, while the absolute agreement is machine
# epsilon).  ``atol`` floors the near-zero shed elements; ``rtol`` bounds the
# O(1) flux (observed relative agreement ~1e-15 — a real divergence would be
# orders of magnitude larger).
_RTOL = 1e-11
_ATOL = 1e-12


def _build_mesh(nx, ny, lvl, ng, bc):
    mesh = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC(bc), bc_xmax=BC(bc), bc_ymin=BC(bc), bc_ymax=BC(bc),
    )
    return SNMesh(
        mesh, Quadrature.level_symmetric(lvl), {0: get_mixture("A", f"{ng}g")}
    )


def _random_inputs(rng, sn_mesh, ng, nx, ny):
    """Het Σ_t + per-ordinate (anisotropic) source + random non-zero inflow."""
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx, ny))
    Q = rng.uniform(0.0, 2.0, size=(sn_mesh.quad.N, ng, nx, ny))
    bf = BoundaryFlux.zeros_on(sn_mesh)
    for face in bf.layout.faces:
        fv = bf.face_view(face)
        fv[...] = rng.uniform(0.0, 1.0, size=fv.shape)
    return sig_t, Q, bf


@pytest.mark.parametrize("nx,ny,lvl,ng,bc", CASES)
def test_scanmarch_sweep_equals_oracle(nx, ny, lvl, ng, bc):
    r"""``ScanMarch.sweep`` ≡ ``FullFieldWavefront.sweep`` to FP-association (nulp).

    Angular flux, scalar flux, AND the post-sweep boundary trace — the
    scan-march reproduces the oracle's lower-triangular solve via the row-march
    schedule (issue #222 G2).  The NON-SQUARE cases catch an x↔y axis swap
    (Mode-2): a square mesh would hide a scan/march-axis confusion.
    """
    rng = np.random.default_rng([nx, ny, lvl, ng])   # deterministic across runs
    sn_mesh = _build_mesh(nx, ny, lvl, ng, bc)
    sig_t, Q, bf_sm = _random_inputs(rng, sn_mesh, ng, nx, ny)
    bf_or = BoundaryFlux.from_mesh(bf_sm.values.copy(), sn_mesh)

    ang_sm, scal_sm = ScanMarch(sn_mesh).sweep(Q, sig_t, bf_sm)
    ang_or, scal_or = FullFieldWavefront(sn_mesh).sweep(Q, sig_t, bf_or)

    np.testing.assert_allclose(ang_sm, ang_or, rtol=_RTOL, atol=_ATOL,
                               err_msg="angular flux")
    np.testing.assert_allclose(scal_sm, scal_or, rtol=_RTOL, atol=_ATOL,
                               err_msg="scalar flux")
    np.testing.assert_allclose(bf_sm.values, bf_or.values, rtol=_RTOL, atol=_ATOL,
                               err_msg="post-sweep boundary trace")


@pytest.mark.parametrize("nx,ny,lvl,ng,bc", [c for c in CASES if c[4] != "vacuum"])
def test_scanmarch_moment_equals_window(nx, ny, lvl, ng, bc):
    r"""ScanMarch moment-OUTPUT ≡ ``MovingFrontierWindow`` moment-output (nulp).

    The Phase-5c moment OUTPUT path (the per-row harmonic projection) reproduces
    the window's per-anti-diagonal projection to FP-association.  The full
    angular field is never materialized; the second return slot is ``None`` (the
    scalar IS ``φ_0^0 = moments[0, 0]``).  ``ℓ≥1`` must carry signal
    (anti-degeneracy) — else the projection-order equivalence is vacuous.  The
    window is a verified moment producer (S4/5c), so this transitively pins the
    scan-march moment output.  ``pytest.fail`` (not bare ``assert``) keeps the
    guards live under ``python -O`` (Mode 8).
    """
    Lm = 1   # P1 — ℓ≥1 load-bearing
    rng = np.random.default_rng([nx, ny, lvl, ng, 7])   # deterministic
    sn_mesh = _build_mesh(nx, ny, lvl, ng, bc)
    sig_t, Q, bf_sm = _random_inputs(rng, sn_mesh, ng, nx, ny)
    bf_win = BoundaryFlux.from_mesh(bf_sm.values.copy(), sn_mesh)
    mp = MomentProjection(
        weights=sn_mesh.quad.weights,
        Y=sn_mesh.quad.spherical_harmonics(Lm),
        L=Lm,
    )

    mom_sm, second_sm = ScanMarch(sn_mesh).sweep(
        Q, sig_t, bf_sm, moment_projection=mp,
    )
    mom_win, _ = MovingFrontierWindow(sn_mesh).sweep(
        Q, sig_t, bf_win, moment_projection=mp,
    )

    if second_sm is not None:
        pytest.fail(
            "ScanMarch moment mode must return None as the 2nd slot (the scalar "
            f"IS φ_0^0 = moments[0,0]); got {type(second_sm).__name__}"
        )
    l0 = float(np.abs(mom_win[0, 0]).max())
    l_ge1 = float(np.abs(mom_win[1:]).max())
    if not l_ge1 > 1e-3 * l0:
        pytest.fail(
            f"ℓ≥1 content degenerate ({l_ge1:.2e} vs ℓ0 {l0:.2e}) — the moment "
            "projection-order equivalence is inert"
        )

    np.testing.assert_allclose(mom_sm, mom_win, rtol=_RTOL, atol=_ATOL,
                               err_msg="moment tensor")
    np.testing.assert_allclose(bf_sm.values, bf_win.values, rtol=_RTOL, atol=_ATOL,
                               err_msg="post-sweep boundary trace")


@pytest.mark.parametrize("nx,ny,lvl,ng,bc", CASES)
def test_scanmarch_residual_equals_oracle(nx, ny, lvl, ng, bc):
    r"""``ScanMarch.residual`` ≡ ``FullFieldWavefront.residual`` (the matvec, G2.c).

    L21 — matvec and sweep are different applications of the SAME operator.  The
    row-march APPLY path (:meth:`ScanMarch.loss_action`, the matvec walk that
    since S6.3 lives on the loss representation, off the operator)
    reconstructs the interior faces from the probe ψ̄ via the ``α = −1``
    reflection scan and evaluates the SAME ``(L+C)ψ`` residual the anti-diagonal
    oracle does — pinned on the **bulk** residual AND the **O.4b boundary-block**
    residual (the latter is what catches a face-shed convention drift).
    Principled-equivalent (the row-march reconstructs the faces in a different
    order), so ``assert_allclose`` (not ``array_equal``) — see ``_RTOL``/``_ATOL``.
    """
    rng = np.random.default_rng([nx, ny, lvl, ng, 13])   # deterministic
    sn_mesh = _build_mesh(nx, ny, lvl, ng, bc)
    sig_t = rng.uniform(0.3, 3.0, size=(ng, nx, ny))
    L = StreamingOperator(sn_mesh, sig_t)

    state = TimedFullField.zeros(
        bulk=AngularFlux, boundary=BoundaryFlux, mesh=sn_mesh,
    )
    state.bulk.values[...] = rng.uniform(-1.0, 1.0, size=state.bulk.values.shape)
    for face in state.boundary.layout.faces:
        fv = state.boundary.face_view(face)
        fv[...] = rng.uniform(0.0, 1.0, size=fv.shape)

    out_sm = ScanMarch(sn_mesh).loss_action(L, state)
    out_or = FullFieldWavefront(sn_mesh).loss_action(L, state)

    np.testing.assert_allclose(
        out_sm.bulk.values, out_or.bulk.values, rtol=_RTOL, atol=_ATOL,
        err_msg="bulk residual",
    )
    np.testing.assert_allclose(
        out_sm.boundary.values, out_or.boundary.values, rtol=_RTOL, atol=_ATOL,
        err_msg="boundary-block residual",
    )
