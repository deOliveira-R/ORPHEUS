r"""L1 MMS tests for SN — 2D Cartesian (1-group and 2-group).

Phase 3.1 (1-group) and Phase 3.2 (2-group heterogeneous) of the
verification campaign.  Both tests verify the **spatial** convergence
of the 2D wavefront diamond-difference sweep
(:func:`orpheus.sn.sweep._sweep_2d_wavefront`) on fixed-source
problems with closed-form reference fluxes.

The 1-group test is the simpler baseline; the 2-group test adds
smooth :math:`\Sigma(x, y)` cross sections with downscatter coupling
to exercise the multigroup scatter assembly in 2D.

See :doc:`/theory/discrete_ordinates` (2D Cartesian MMS sections)
for the full derivation and :mod:`orpheus.derivations.sn_mms` for
the reference solution.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations import continuous_get
from orpheus.derivations.continuous.mms.sn import (
    build_2d_cartesian_mms_case,
    build_2d_cartesian_heterogeneous_mms_case,
)
from orpheus.sn import solve_sn_fixed_source


def _l2_2d(err: np.ndarray, volumes: np.ndarray) -> float:
    r"""Volume-weighted :math:`L^{2}` norm on a 2D mesh."""
    return float(np.sqrt(np.sum(volumes * err * err)))


# ═══════════════════════════════════════════════════════════════════════
# Phase 3.1 — 1-group 2D Cartesian MMS
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l1
@pytest.mark.verifies(
    "transport-cartesian-2d", "dd-cartesian-2d",
    "sn-mms-2d-psi", "sn-mms-2d-qext",
)
def test_sn_2d_cartesian_mms_converges_second_order():
    r"""Diamond-difference SN on a 2D Cartesian mesh shows :math:`\mathcal{O}(h^2)`.

    Runs the separable MMS problem on four refinements of a square
    mesh (nx = ny = 10, 20, 40, 80) and asserts the observed
    convergence order between successive refinements is > 1.9.
    """
    case = build_2d_cartesian_mms_case()

    n_cells = [10, 20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            boundary_condition="vacuum",
            max_inner=500,
            inner_tol=1e-13,
        )
        phi_num = result.scalar_flux[:, :, 0]         # (nx, ny)
        phi_ref = case.phi_exact(mesh.centers_x, mesh.centers_y)
        errors.append(_l2_2d(phi_num - phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2) convergence, got orders={orders} "
        f"from errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3


# ═══════════════════════════════════════════════════════════════════════
# Phase 3.2 — 2-group heterogeneous 2D Cartesian MMS
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.l1
@pytest.mark.verifies(
    "transport-cartesian-2d", "dd-cartesian-2d",
    "multigroup", "mg-balance",
    "sn-mms-2d-2g-psi", "sn-mms-2d-2g-qext",
)
def test_sn_2d_cartesian_2g_heterogeneous_mms_converges_second_order():
    r"""2-group heterogeneous SN on a 2D Cartesian mesh shows :math:`\mathcal{O}(h^2)`.

    Smooth :math:`\Sigma(x, y)` cross sections with downscatter coupling,
    verified at four mesh refinements.  Both groups must converge at
    :math:`> 1.9` independently.
    """
    case = build_2d_cartesian_heterogeneous_mms_case()

    n_cells = [10, 20, 40, 80]
    errs_g0 = []
    errs_g1 = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        materials = case.build_materials(mesh)
        Q = case.external_source(mesh)

        result = solve_sn_fixed_source(
            materials, mesh, case.quadrature, Q,
            boundary_condition="vacuum",
            max_inner=500, inner_tol=1e-12,
        )

        phi_solver = result.scalar_flux  # (nx, ny, ng)
        for g, errs in [(0, errs_g0), (1, errs_g1)]:
            phi_ref = case.phi_exact(mesh.centers_x, mesh.centers_y, g)
            errs.append(_l2_2d(phi_solver[:, :, g] - phi_ref, mesh.volumes))

    errs_g0_arr = np.asarray(errs_g0)
    errs_g1_arr = np.asarray(errs_g1)
    orders_g0 = np.log2(errs_g0_arr[:-1] / errs_g0_arr[1:])
    orders_g1 = np.log2(errs_g1_arr[:-1] / errs_g1_arr[1:])

    assert np.all(orders_g0 > 1.9), (
        f"Group 0: errors={errs_g0_arr}, orders={orders_g0}"
    )
    assert np.all(orders_g1 > 1.9), (
        f"Group 1: errors={errs_g1_arr}, orders={orders_g1}"
    )
    assert 1e-7 < errs_g0_arr[-1] < 1e-3
    assert 1e-7 < errs_g1_arr[-1] < 1e-3


@pytest.mark.l1
@pytest.mark.verifies("sn-mms-2d-2g-psi")
def test_sn_2d_heterogeneous_mms_positive_absorption_everywhere():
    r"""The 2D smooth Σ(x,y) must give Σ_a > 0 everywhere.

    Regression guard on the 2D cross-section functions: any tweak
    that makes absorption go negative would produce unphysical
    materials.
    """
    case = build_2d_cartesian_heterogeneous_mms_case()
    Lx, Ly = case.length_x, case.length_y

    x = np.linspace(0.0, Lx, 101)
    y = np.linspace(0.0, Ly, 101)
    xx, yy = np.meshgrid(x, y, indexing="ij")
    xf, yf = xx.ravel(), yy.ravel()

    for g in range(2):
        sig_t = case.sigma_t_fn(xf, yf, g)
        sig_s_out = sum(
            case.sigma_s_fn(xf, yf, g, gp) for gp in range(2)
        )
        sig_a = sig_t - sig_s_out
        assert np.all(sig_a > 0), (
            f"Σ_a,{g}(x,y) non-positive: min = {np.min(sig_a):.4e}"
        )
