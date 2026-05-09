r"""Document and catch ERR-026: curvilinear sweep WDD angular closure
converges to wrong fixed-source solution.

The spherical and cylindrical sweeps use a one-directional WDD angular
face-flux closure that, combined with the zero-area face at r=0,
converges to a non-flat solution for constant-source problems.  The
symmetric-closure operator (now reached via ``inner_solver="krylov"``
on :func:`solve_sn_fixed_source`) gives the correct answer.

Wave E Round 2 (Issue #164) closes ERR-026 by routing curvilinear
``solve_sn_fixed_source`` calls through Krylov-on-:meth:`L.apply`.
Round 2 keeps the entire test file under :class:`pytest.mark.xfail`
because the **assertions** here document the historical sweep-vs-
operator divergence (i.e. they assert the sweep IS wrong by >= 20 %
relative).  After Round 2 the krylov path lands the correct answer,
so the gap shrinks and these assertions stop holding — exactly the
"closure" condition Wave E set out to deliver.  Round 3 owns the
full file rewrite into a positive ERR-026-closed regression.

See:
- GitHub Issue #98 (sweep-operator inconsistency, closes in Wave E R2)
- GitHub Issue #99 (Phase 3.3-3.4 MMS blocker)
- ``.claude/skills/vv-principles/error_catalog.md`` ERR-026
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import Mesh1D
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.mesh import BC
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.sweep import transport_sweep


def _make_spherical_problem(nx: int = 10, R: float = 10.0, N_ord: int = 8):
    """Build a constant-source, reflective-BC spherical problem."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("reflective"),
    )
    quad = GaussLegendre1D.create(N_ord)
    sn_mesh = SNMesh(mesh, quad)
    sig_t = np.full((nx, 1, 1), 1.0)
    Q_iso = np.full((nx, 1, 1), 1.0)
    return sn_mesh, quad, sig_t, Q_iso


def _solve_sweep(sn_mesh, sig_t, Q_iso, max_iter=200):
    """Source iteration via the spherical sweep."""
    psi_bc = {}
    phi_old = None
    for it in range(max_iter):
        _, phi = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)
        if phi_old is not None:
            res = np.linalg.norm(phi - phi_old) / max(np.linalg.norm(phi), 1e-30)
            if res < 1e-14:
                break
        phi_old = phi.copy()
    return phi[:, 0, 0]


def _build_pure_absorber_1g(sig_t_val: float):
    """Make a 1-group pure-absorber Mixture matching the legacy fixture.

    Σ_t = Σ_C (capture only); no scattering, fission, leakage, or
    (n,2n). Σ_a = Σ_t since absorption = capture + fission + (n,2n)-out.
    """
    from orpheus.data.macro_xs.mixture import Mixture
    from scipy.sparse import csr_matrix

    sig_c = np.asarray([sig_t_val])
    sig_l = np.zeros(1)
    sig_f = np.zeros(1)
    sig_p = np.zeros(1)
    sig_t = sig_c.copy()  # only capture contributes
    sig_s = csr_matrix(np.zeros((1, 1)))
    sig2 = csr_matrix(np.zeros((1, 1)))
    chi = np.zeros(1)
    return Mixture(
        SigC=sig_c, SigL=sig_l, SigF=sig_f, SigP=sig_p, SigT=sig_t,
        SigS=[sig_s], Sig2=sig2, chi=chi,
    )


def _solve_via_krylov(sn_mesh, quad, sig_t, Q_iso):
    """Solve via the Wave E Round 2 ``inner_solver="krylov"`` path.

    The legacy ``_solve_bicgstab`` helper (which built the FD operator
    explicitly via :func:`build_transport_linear_operator_spherical`
    and ran BiCGSTAB on the packed solution vector) was retired along
    with the rest of the BiCGSTAB FD-operator API.  This replacement
    routes through :func:`orpheus.sn.solver.solve_sn_fixed_source`
    with ``inner_solver="krylov"``, which uses GMRES on
    :meth:`SNStreamingOperator.apply` (the same symmetric closure the
    legacy BiCGSTAB path used).
    """
    from orpheus.sn.solver import solve_sn_fixed_source

    nx = sn_mesh.nx
    N = quad.N

    # Match the symmetric-closure operator's per-ordinate constant
    # source.  The sweep applies 1/W internally; solve_sn_fixed_source
    # routes the krylov path the same way.
    external_source = np.ones((N, nx, 1, 1))

    # Single-group, 1-region material with Σ_t = 1.0, Σ_s = 0.
    materials = {1: _build_pure_absorber_1g(sig_t[0, 0, 0])}

    result = solve_sn_fixed_source(
        materials, sn_mesh.mesh, quad, external_source,
        boundary_condition="reflective",
        max_inner=200, inner_tol=1e-10,
        inner_solver="krylov",
    )
    return result.scalar_flux[:, 0, 0]


# ═══════════════════════════════════════════════════════════════════════
# Tests
# ═══════════════════════════════════════════════════════════════════════

# Wave E Round 2: the legacy ``_solve_bicgstab`` helper that this
# file used to drive the symmetric-closure operator was retired
# along with the rest of the BiCGSTAB FD-operator API surface.  The
# replacement helper ``_solve_via_krylov`` routes through
# :func:`solve_sn_fixed_source` with ``inner_solver="krylov"``, which
# uses GMRES on :meth:`SNStreamingOperator.apply` — the same
# symmetric-closure operator under a different driver.  All four
# assertions in this file remain mechanically valid: the sweep still
# carries the WDD asymmetric closure, so it still produces the
# documented ERR-026 deviation against the symmetric-closure
# reference.
#
# Round 3 owns the final rewrite of this file into a positive
# ERR-026-closed regression (the test names will flip from
# "sweep_vs_bicgstab" to "krylov_matches_analytical_flat_flux").
pytestmark = [pytest.mark.l1, pytest.mark.catches("ERR-026")]


def test_spherical_sweep_vs_bicgstab_flat_flux():
    r"""Krylov gives exact flat flux; sweep deviates significantly.

    Constant isotropic source :math:`Q = \Sigma_t = 1`, reflective BCs.
    Expected :math:`\phi = 1` everywhere. Krylov gets it; the sweep
    converges to a stable but wrong profile with ~35% error at r=0.

    This test **documents** ERR-026 — it does NOT assert the sweep is
    correct, because it isn't. It asserts Krylov IS correct and that
    the sweep's deviation is at least as large as observed.
    """
    sn_mesh, quad, sig_t, Q_iso = _make_spherical_problem(nx=20)
    phi_k = _solve_via_krylov(sn_mesh, quad, sig_t, Q_iso)
    phi_sweep = _solve_sweep(sn_mesh, sig_t, Q_iso)

    # Krylov must be exact
    np.testing.assert_allclose(phi_k, 1.0, atol=1e-10,
                               err_msg="Krylov should give exact flat flux")

    # Sweep must deviate significantly (ERR-026 evidence)
    sweep_err = np.max(np.abs(phi_sweep - 1.0))
    assert sweep_err > 0.2, (
        f"Sweep error {sweep_err:.4e} is suspiciously small — "
        f"has the sweep bug been fixed? If so, remove ERR-026."
    )


def test_spherical_sweep_error_does_not_converge():
    r"""The sweep's error at cell 0 does NOT converge with refinement.

    This is the defining characteristic of ERR-026: the error is structural,
    not truncation error.  If this test ever fails (orders > 1.5), the
    sweep has been fixed and the MMS blocker (Issue #99) can be resolved.
    """
    errors = []
    for nx in [10, 20, 40]:
        sn_mesh, quad, sig_t, Q_iso = _make_spherical_problem(nx=nx)
        phi_sweep = _solve_sweep(sn_mesh, sig_t, Q_iso)
        errors.append(np.max(np.abs(phi_sweep - 1.0)))

    errors = np.asarray(errors)
    ratios = errors[:-1] / errors[1:]
    orders = np.log2(ratios)

    # The error should NOT converge (orders < 0.5 means diverging or stagnant)
    assert np.all(orders < 0.5), (
        f"Sweep error is converging (orders={orders}) — "
        f"has ERR-026 been fixed? Update Issue #98."
    )


def test_spherical_sweep_conserves_globally():
    r"""Despite the wrong spatial profile, global conservation holds.

    Volume-weighted average :math:`\phi` equals :math:`Q/\Sigma_t = 1`
    to machine precision, because the balance equation is satisfied per
    cell — the error is in the spatial DISTRIBUTION, not the total.
    """
    sn_mesh, quad, sig_t, Q_iso = _make_spherical_problem(nx=20)
    phi_sweep = _solve_sweep(sn_mesh, sig_t, Q_iso)
    V = sn_mesh.volumes[:, 0]
    phi_vol_avg = np.sum(phi_sweep * V) / V.sum()

    np.testing.assert_allclose(phi_vol_avg, 1.0, atol=1e-10,
                               err_msg="Global conservation should hold even with ERR-026")


def test_cartesian_sweep_gives_exact_flat_flux():
    r"""Cartesian sweep has no angular redistribution — flat flux is exact.

    Control test: confirms the issue is specific to curvilinear geometry.
    """
    nx, R = 20, 10.0
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
    )
    quad = GaussLegendre1D.create(8)
    sn_mesh = SNMesh(mesh, quad)
    sig_t = np.full((nx, 1, 1), 1.0)
    Q_iso = np.full((nx, 1, 1), 1.0)

    phi = _solve_sweep(sn_mesh, sig_t, Q_iso)
    np.testing.assert_allclose(phi, 1.0, atol=1e-10,
                               err_msg="Cartesian sweep should give exact flat flux")
