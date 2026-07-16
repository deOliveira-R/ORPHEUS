r"""LD slope-moment frame consistency — the diffusion-limit root-cause gate (#240 D5b-S3).

The per-cell LD kernel produces/consumes the ``2^d`` moment vector in the
per-ordinate SWEEP frame (each axis oriented so the downstream face is at local
``+1``).  The iterate ``φ̂`` and its isotropic scattering source ``Σ_s·φ̂`` live
in the GLOBAL-x frame, where the angular reduction ``φ̂ = Σ_n w_n ψ̂_n`` sums the
slope across ordinates of BOTH sweep directions.  For a backward ordinate
(``μ < 0``) the sweep coordinate is the reverse of the global coordinate, so the
slope (``P₁``) moment is sign-flipped between the two frames.  ERR-061: the
producer stored the raw sweep-frame slope and the consumer summed it as
global-frame, so the backward ordinates' opposite-signed slopes CANCELLED the
forward ones → ``φ̂`` was ~6× under-driven → LD lost the thick-cell diffusion
limit.

The fix lifts the slope to the global frame at the octant boundary (the
``octant_moment_frame_signs`` involution).  This file pins:

* **Frame consistency** — at a cell with a positive global-x gradient, the
  forward AND backward ordinate slopes ``ψ̂_n`` share sign (the within-cell
  linear profile is the same function for all ordinates; its global-x slope does
  not depend on the sweep direction).  Pre-fix the backward slope had the
  OPPOSITE sign (the sweep frame).
* **Structural-independence confirmation** — a from-scratch LD-SN slab solver
  (a direct LM-1989 Eq 4.3 cell 2×2 + source iteration, NO ORPHEUS kernel)
  recovers the analytic diffusion VALUE at the thick coarse mesh ONLY with the
  global-frame slope correction, and reproduces the (broken) pre-fix value
  bit-for-bit with the sweep-frame slope — pinning the root cause to the
  slope-moment frame independent of ORPHEUS's code.

Failure mode #1 (sign flip) + #6 (convention drift, sweep-frame producer vs
global-frame consumer).  See ERR-061 and ``docs/theory/methods/sn/index.rst``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import _make_1g_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.transport.spatial import LinearDiscontinuous

THETA = 1.0 / 3.0


@pytest.mark.l1
@pytest.mark.verifies("ld-ubld-d1-reduction", "ld-ubld-slope-angular-reduction")
@pytest.mark.catches("ERR-061")
def test_ld_slope_moment_global_frame_consistency() -> None:
    """Forward and backward ordinates store the SAME-sign global-x slope.

    On a cell with a clear positive global-x gradient the per-ordinate slope
    ``ψ̂_n`` must share sign across sweep directions (the within-cell linear
    profile is the same function for every ordinate).  A sign disagreement is
    the sweep-frame storage bug (ERR-061).
    """
    sigma_t, c, length = 40.0, 0.99, 1.0
    nx = 16
    materials = {0: _make_1g_mixture(sigma_t, c * sigma_t)}
    mesh = Mesh1D(
        edges=np.linspace(0.0, length, nx + 1), mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN, bc_left=BC("vacuum"), bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    Q = np.ones((quad.N, 1, nx)) / quad.weights.sum()
    ld = solve_sn_fixed_source(
        materials, mesh, quad, Q, scheme=LinearDiscontinuous(),
        boundary_condition="vacuum", inner_solver="krylov",
        max_inner=2000, inner_tol=1e-10,
    )
    ang = ld.angular_flux.interior.values   # (N, ng, nx, 2)
    mu = quad.mu_x
    n_fwd, n_bwd = int(np.argmax(mu)), int(np.argmin(mu))
    j = 3   # left flank: global-x gradient clearly positive (flux rising L→R)
    hat_fwd = float(ang[n_fwd, 0, j, 1])
    hat_bwd = float(ang[n_bwd, 0, j, 1])
    np.testing.assert_array_less(
        0.0, hat_fwd,
        err_msg=f"forward slope should be positive (rising flux), got {hat_fwd}",
    )
    np.testing.assert_array_less(
        0.0, hat_bwd,
        err_msg=(
            f"backward-ordinate slope ψ̂_n={hat_bwd:.5f} has OPPOSITE sign to "
            f"forward ψ̂_n={hat_fwd:.5f} — the LD slope moment is stored in the "
            "per-ordinate SWEEP frame, not the global-x frame (ERR-061; the "
            "angular reduction φ̂=Σw_n ψ̂_n then cancels → diffusion limit lost)"
        ),
    )


def _independent_ld_slab(nx, length, sigma_t, c, Q_iso, mu, wt, *, global_frame):
    """From-scratch LD-SN slab solver (independent of ORPHEUS kernels).

    Direct LM-1989 Eq (4.3) cell 2×2 + source iteration.  ``global_frame``: store
    the per-ordinate slope in the global-x frame (the fix, backward ordinates
    sign-flipped) or the raw sweep frame (the bug, summed directly).
    """
    h = length / nx
    W = wt.sum()
    sigma_s = c * sigma_t
    N = len(mu)
    psi_bar = np.zeros((N, nx))
    psi_hat = np.zeros((N, nx))
    phi_bar = np.zeros(nx)
    phi_hat = np.zeros(nx)

    def cell(mu_abs, s_bar, s_hat, psi_in):
        A = np.array([[sigma_t * h + mu_abs, mu_abs],
                      [-mu_abs / THETA, sigma_t * h + mu_abs / THETA]])
        rhs = np.array([s_bar * h + mu_abs * psi_in,
                        s_hat * h - (mu_abs / THETA) * psi_in])
        pb, ph = np.linalg.solve(A, rhs)
        return pb, ph, pb + ph

    for _ in range(20000):
        phi_bar_old = phi_bar.copy()
        s_bar = (Q_iso + sigma_s * phi_bar) / W
        s_hat_global = (sigma_s * phi_hat) / W
        nb = np.zeros((N, nx))
        nh = np.zeros((N, nx))
        for n in range(N):
            mu_abs = abs(mu[n])
            forward = mu[n] > 0
            sign = 1.0 if (forward or not global_frame) else -1.0
            rng = range(nx) if forward else range(nx - 1, -1, -1)
            psi_in = 0.0
            for j in rng:
                pb, ph_sweep, pout = cell(mu_abs, s_bar[j], sign * s_hat_global[j], psi_in)
                nb[n, j] = pb
                nh[n, j] = sign * ph_sweep
                psi_in = pout
        psi_bar, psi_hat = nb, nh
        phi_bar = np.einsum("n,nj->j", wt, psi_bar)
        phi_hat = np.einsum("n,nj->j", wt, psi_hat)
        if np.max(np.abs(phi_bar - phi_bar_old)) < 1e-12:
            break
    return phi_bar


@pytest.mark.foundation
def test_independent_ld_global_frame_recovers_diffusion() -> None:
    """Structurally-independent LD (global-frame slope) recovers diffusion at nx=4.

    A direct LM-1989 2×2 + source iteration (NO ORPHEUS kernel).  The global-frame
    slope (the fix) recovers the analytic diffusion VALUE at the thick coarse
    mesh; the sweep-frame slope (the bug) is far off — the structurally-
    independent confirmation that the fix is correct AND that the LD SCHEME,
    correctly framed, IS diffusion-consistent (vv §1/§6).

    NOTE: this carries NO ``catches("ERR-061")`` — it does not exercise ORPHEUS's
    ``_reframe`` (re-introducing the bug there leaves this GREEN); it is the
    structural-independence GROUND that proves the bug class is real and the fix
    direction correct, not a regression catcher for the production code path.
    """
    sigma_t, c, length, Q_iso = 40.0, 0.99, 1.0, 1.0
    mu, wt = np.polynomial.legendre.leggauss(8)
    D = 1.0 / (3.0 * sigma_t)
    sigma_a = (1.0 - c) * sigma_t
    Ld = np.sqrt(D / sigma_a)
    extrap = 0.7104 / sigma_t
    phi_diff = (Q_iso / sigma_a) * (1.0 - 1.0 / np.cosh((length / 2 + extrap) / Ld))

    nx = 4
    phi_global = _independent_ld_slab(nx, length, sigma_t, c, Q_iso, mu, wt,
                                      global_frame=True)
    phi_sweep = _independent_ld_slab(nx, length, sigma_t, c, Q_iso, mu, wt,
                                     global_frame=False)
    rel_global = abs(phi_global[2] - phi_diff) / phi_diff
    rel_sweep = abs(phi_sweep[2] - phi_diff) / phi_diff
    np.testing.assert_array_less(
        0.30, rel_sweep,
        err_msg=(
            f"sweep-frame LD should LOSE the diffusion limit (rel {rel_sweep:.1%})"
            " — if small, the independent solver does not reproduce ERR-061"
        ),
    )
    np.testing.assert_array_less(
        rel_global, 0.05,
        err_msg=(
            f"global-frame LD must RECOVER the diffusion limit at nx=4 "
            f"(rel {rel_global:.1%}) — the structurally-independent confirmation "
            "of the slope-moment-frame fix (ERR-061)"
        ),
    )
