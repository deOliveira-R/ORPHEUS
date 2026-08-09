"""Diagnostic: source-iteration cost on an ALL-reflective pure absorber
scales as ~1/Sigma_t and explodes with DIMENSION -- so the entry's fixed
``max_inner=1000`` is a dimension-blind cap.

Created by numerics-investigator on 2026-08-08.

The general property behind the
``test_d3_pure_absorber_per_ordinate_psi_exact`` red.  With NO scattering and
NO leakage, the only coupling left is the reflective boundary and the only
damping is absorption.  The slow mode is the DD face sawtooth: it has ZERO
cell average by construction, so the collision term ``Sigma_t V psi_c``
barely sees it, and it decays only through the inter-axis balance mismatch.

`[M]` 2026-08-08, ``inner_tol=1e-13``, level-symmetric S4, boundary-G-S:

    d=1 slab (3 cells)          n_inner =   32
    d=2 box   (3,4)             n_inner =  258
    d=3 box   (3,4,5)           n_inner = 1631      <-- past max_inner=1000

    d=3, Sigma_t = 0.4 / 0.8    n_inner = 3093
    d=3, Sigma_t = 0.8 / 1.6    n_inner = 1631
    d=3, Sigma_t = 1.6 / 3.2    n_inner =  850
    d=3, Sigma_t = 3.2 / 6.4    n_inner =  437
    (Sigma_t * n_inner ~ 1240..1400 -- absorption-limited, as the mechanism
     predicts)

These are the two gates a permanent home should keep; they are properties of
the METHOD, not of one fixture, and they are what makes the fixed cap wrong.

If this catches a real bug, promote to ``tests/sn/solve/`` (a new
``test_reflective_si_iteration_budget.py``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.mesh.axis import AxisMesh

_Q_G = np.array([1.0, 0.5])
_BUDGET = 8000          # generous headroom so every case truly converges


def _absorber(sig):
    sig = np.asarray(sig, float)
    return make_mixture(
        sig_t=sig, sig_c=sig, sig_f=np.zeros(2), nu=np.zeros(2),
        chi=np.zeros(2), sig_s=np.zeros((2, 2)),
    )


def _run(mesh, quad, spatial_shape, sig):
    mix = _absorber(sig)
    W = float(np.sum(quad.weights))
    q = np.broadcast_to(
        (_Q_G / W).reshape(1, 2, *([1] * len(spatial_shape))),
        (quad.N, 2, *spatial_shape),
    ).copy()
    sol = solve_sn_fixed_source(
        {0: mix}, mesh, quad, external_source=q,
        boundary_condition="reflective", inner_tol=1e-13,
        max_inner=_BUDGET,
    )
    psi = np.asarray(sol.angular_flux.interior.values)
    err = max(
        float(np.max(np.abs(psi[:, g] - _Q_G[g] / (W * sig[g]))
                     / (_Q_G[g] / (W * sig[g]))))
        for g in range(2)
    )
    return sol.history, err


def _axes3():
    return tuple(
        AxisMesh(edges=np.linspace(0.0, ext, n + 1), bc_low=None, bc_high=None)
        for ext, n in zip((1.0, 2.0, 3.0), (3, 4, 5))
    )


@pytest.mark.l1
@pytest.mark.slow
def test_reflective_absorber_iteration_count_grows_with_dimension() -> None:
    """d=1 -> d=2 -> d=3 costs ~32 -> ~258 -> ~1631 sweeps at 1e-13.

    The claim under test is ORDERING, not the exact counts (those move with
    the quadrature and the mesh): each added dimension multiplies the
    reflective-SI budget by roughly an order, so d=3 crosses the entry's
    ``max_inner=1000`` while d=1/d=2 sit comfortably below it.  That is
    precisely why a fixed cap cannot serve all three.
    """
    sig = [0.8, 1.6]
    quad = Quadrature.level_symmetric(sn_order=4)

    h1, e1 = _run(
        Mesh1D(edges=np.linspace(0.0, 1.0, 4),
               mat_ids=np.zeros(3, dtype=int), coord=CoordSystem.CARTESIAN),
        quad, (3,), sig,
    )
    h2, e2 = _run(
        Mesh2D(edges_x=np.linspace(0.0, 1.0, 4),
               edges_y=np.linspace(0.0, 2.0, 5),
               mat_map=np.zeros((3, 4), dtype=int)),
        quad, (3, 4), sig,
    )
    h3, e3 = _run(_axes3(), quad, (3, 4, 5), sig)

    for tag, h, e in (("d=1", h1, e1), ("d=2", h2, e2), ("d=3", h3, e3)):
        np.testing.assert_equal(bool(h.converged), True)
        np.testing.assert_array_less(
            e, 1e-10, err_msg=f"{tag}: converged answer must be the flat field",
        )

    # the ordering, with margin
    np.testing.assert_array_less(h1.n_inner, h2.n_inner)
    np.testing.assert_array_less(h2.n_inner, h3.n_inner)
    # d=1 and d=2 fit inside the entry default; d=3 does NOT.
    np.testing.assert_array_less(h1.n_inner, 1000)
    np.testing.assert_array_less(h2.n_inner, 1000)
    if h3.n_inner <= 1000:
        pytest.fail(
            f"d=3 now converges in {h3.n_inner} <= 1000 sweeps — the "
            f"iteration budget changed (an accelerator? a better splitting?). "
            f"Re-derive this gate and the max_inner guidance with it."
        )


@pytest.mark.l1
@pytest.mark.slow
def test_reflective_absorber_iteration_count_is_absorption_limited() -> None:
    """``Sigma_t * n_inner`` is ~invariant: the slow mode is damped ONLY by
    absorption.

    This is the mechanism pin.  The DD face sawtooth has zero cell average,
    so it is invisible to every angular/spatial average; what removes it is
    the absorption suffered along the reflective loop.  Doubling Sigma_t
    therefore halves the sweep count.  If this product stops being invariant,
    the damping mechanism changed and every ``max_inner`` heuristic derived
    from it is void.
    """
    quad = Quadrature.level_symmetric(sn_order=4)
    products = []
    for s in (0.8, 1.6, 3.2):
        h, e = _run(_axes3(), quad, (3, 4, 5), [s, 2 * s])
        np.testing.assert_equal(bool(h.converged), True)
        np.testing.assert_array_less(e, 1e-10)
        products.append(s * h.n_inner)

    products = np.asarray(products, float)
    spread = float(products.max() / products.min())
    if spread > 1.35:
        pytest.fail(
            f"Sigma_t * n_inner is not invariant (spread {spread:.3f} > 1.35, "
            f"values {products}) — the reflective-SI slow mode is no longer "
            f"absorption-limited."
        )
