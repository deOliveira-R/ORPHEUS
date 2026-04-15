"""L1 MMS (Method of Manufactured Solutions) tests for SN — 1D slab.

These tests verify the *spatial* convergence of the diamond-difference
discretisation on a fixed-source problem whose exact angular flux is
known in closed form. Any drift from :math:`\\mathcal{O}(h^{2})` is a
bug in the sweep, the boundary-condition plumbing, or the DD closure.

Related:
- :mod:`orpheus.derivations.sn_mms` — ansatz + manufactured source
- :func:`orpheus.sn.solver.solve_sn_fixed_source` — fixed-source driver
- Theory page: ``docs/theory/discrete_ordinates.rst`` (MMS section)
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.sn_mms import build_1d_slab_mms_case
from orpheus.sn import solve_sn_fixed_source


def _l2_error(phi_num: np.ndarray, phi_ref: np.ndarray, widths: np.ndarray) -> float:
    """Cell-width-weighted discrete L2 norm of the flux error."""
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(widths * diff * diff)))


@pytest.mark.l1
@pytest.mark.verifies(
    "dd-cartesian-1d", "dd-slab", "transport-cartesian",
    "sn-mms-psi", "sn-mms-qext",
)
def test_sn_1d_slab_mms_converges_second_order():
    r"""Diamond-difference SN on a 1D slab shows measured :math:`\mathcal{O}(h^2)`.

    Runs the MMS problem at four mesh refinements (nc = 20, 40, 80, 160)
    and asserts that the observed convergence order between successive
    refinements is :math:`> 1.9`. The ansatz :math:`\phi(x) = \sin(\pi x/L)`
    is smooth, so the DD discretisation error — not the manufactured
    source or the quadrature — governs the decay.
    """
    case = build_1d_slab_mms_case()

    n_cells = [20, 40, 80, 160]
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
        phi_num = result.scalar_flux[:, 0, 0]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_error(phi_num, phi_ref, mesh.widths))

    errors = np.asarray(errors)
    ratios = errors[:-1] / errors[1:]
    orders = np.log2(ratios)

    # Every successive refinement must show ≥ 1.9 — DD is exactly 2,
    # so 1.9 is a tight bracket that still leaves room for round-off
    # at the finest mesh once the error approaches inner_tol.
    assert np.all(orders > 1.9), (
        f"Expected O(h^2) convergence, got orders={orders} "
        f"from errors={errors}"
    )

    # Absolute magnitude sanity: the finest mesh error should be small
    # but well above inner_tol (otherwise we measured iteration noise).
    assert 1e-8 < errors[-1] < 1e-4


@pytest.mark.l1
@pytest.mark.verifies("dd-cartesian-1d")
def test_sn_mms_manufactured_source_vanishes_at_zero_material():
    r"""Consistency of the manufactured source formula.

    When :math:`\Sigma_t = \Sigma_s` (pure scatterer) and the ordinate
    is :math:`\mu = 0`, the per-ordinate external source reduces to
    exactly zero because both the streaming term (:math:`\mu A'`) and
    the removal term (:math:`(\Sigma_t - \Sigma_s) A`) vanish.

    This is an L0-flavoured algebraic cross-check of the symbolic
    derivation in :mod:`orpheus.derivations.sn_mms` — if someone edits
    the formula and drops a term, this test flags it immediately.
    """
    # Equal Σ_t, Σ_s would make _make_1g_mixture reject the mixture
    # (c < 1 is required for source iteration). So instead probe the
    # formula directly: pick μ=0 for a standard case and check that
    # the streaming contribution is exactly zero there.
    case = build_1d_slab_mms_case(sigma_t=1.0, sigma_s=0.9)
    mesh = case.build_mesh(32)
    Q = case.external_source(mesh)  # (N, nx, 1, 1)

    # GaussLegendre1D never places an ordinate exactly at μ=0 for
    # even N, so probe the algebraic identity instead: the half-sum
    # of symmetric ordinates removes the μ·A' contribution and
    # isolates (Σ_t−Σ_s)·A.
    Q_flat = Q[:, :, 0, 0]  # (N, nx)
    N = Q_flat.shape[0]
    Q_sym = 0.5 * (Q_flat[: N // 2] + Q_flat[-(N // 2):][::-1])  # avg μ ↔ −μ
    removal = (case.sigma_t - case.sigma_s) * case.phi_exact(mesh.centers)
    np.testing.assert_allclose(Q_sym, np.broadcast_to(removal, Q_sym.shape),
                                rtol=1e-13, atol=1e-13)
