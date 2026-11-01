r"""LD slope-moment frame consistency on the SN solve (ERR-061), moved from
``tests/gates/transport/spatial/test_ld_slope_frame.py`` by layer DAG Phase 4 (2026-11-01).

At a cell with a positive global-x gradient, the forward AND backward ordinate slopes of an SN fixed-source solve
share sign: the within-cell linear profile is one function for all ordinates, so its global-x slope does not depend on
the sweep direction. Pre-fix (ERR-061) the backward slope had the opposite sign (the sweep frame). The row drives
``solve_sn_fixed_source``, an SN entry, so it lives with the SN slab sweep; the transport file keeps the structurally
independent from-scratch LD solver that confirms the root cause, and its module docstring carries the full argument.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import _make_1g_mixture
from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn_fixed_source
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.transport.spatial import LinearDiscontinuous
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry.
_CONDITIONS = {0: WallCondition.vacuum()}


@pytest.mark.rests_on("tests/gates/transport/spatial/test_ubld_reads_the_interval_basis.py::test_the_octant_frame_map_is_the_reflection_signs")
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
    mesh = Mesher(StructuredGeometry.from_homogeneous(length, 0)).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    quad = Quadrature.gauss_legendre(8)
    Q = np.ones((quad.n_ordinates, 1, nx)) / quad.weights.sum()
    ld = solve_sn_fixed_source(
        materials, mesh, quad, Q, scheme=LinearDiscontinuous(),
        inner_solver="krylov",
        max_inner=2000, inner_tol=1e-10, boundary_conditions=conditions_for(mesh, _CONDITIONS),
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
