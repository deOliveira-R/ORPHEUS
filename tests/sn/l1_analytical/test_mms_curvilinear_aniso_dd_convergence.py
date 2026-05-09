r"""L1 MMS convergence tests for curvilinear SN with **anisotropic** ansatz.

Companion to :file:`tests/sn/test_mms_curvilinear.py` (legacy isotropic
ansatz, currently failing on main with order ≈ 0 — ERR-026). This file
exercises the angular-redistribution operator directly via the
:math:`(A(r) + B(r) \cdot \zeta) / W` ansatz from the method-implementer's
Phase-1 work, where :math:`\zeta = \mu` for sphere and :math:`\zeta = \eta`
for cylinder.

Failure-mode #7 from :ref:`vv-principles` motivates the redesign: the
isotropic ansatz nulls the angular redistribution by construction, so
the MMS test cannot detect ERR-026-class bugs (curvilinear sweep WDD
wrong fixed point). The anisotropic case adds the :math:`B(r) \zeta`
term that activates :math:`(1-\mu^2) B/r` (sphere) or :math:`\xi^2 B/r`
(cylinder) — the angular-redistribution analytic — making the
curvilinear sweep's hardest math observable under refinement.

**Both tests are currently expected to fail on main** with the same
order ≈ 0 signature as the legacy isotropic ones (ERR-026 reaches both
ansatz). They are tagged ``@pytest.mark.xfail(strict=True, ...)``: the
SN reshape Issues 11 (``SNStreamingOperator``) and 12 (unified sweep)
must close ERR-026; when they do, these tests pass and ``strict=True``
flips xpass to a failure that forces removal of the xfail marker.

Pairing rationale (the test-architect's failure-narrowing instrument):
the anisotropic and isotropic cases differ ONLY in the ``B(r) \zeta``
term. If both fail with the same order signature, the bug is upstream
of angular redistribution (e.g., DD spatial closure on curvilinear
mesh). If the isotropic passes and the anisotropic fails, the bug is
in the angular-redistribution code path. If both pass, the curvilinear
sweep is sound. The pair is the diagnostic.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_spherical_anisotropic_mms_case,
    build_cylindrical_anisotropic_mms_case,
)
from orpheus.sn import solve_sn_fixed_source


pytestmark = pytest.mark.l1


def _l2_1d(phi_num: np.ndarray, phi_ref: np.ndarray, volumes: np.ndarray) -> float:
    """Volume-weighted L2 norm for 1D curvilinear meshes."""
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(volumes * diff * diff)))


@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-026 — curvilinear sweep WDD wrong fixed point. The legacy "
        "isotropic spherical MMS test in tests/sn/test_mms_curvilinear.py "
        "currently fails on main with order ≈ 0 instead of O(h²); the "
        "anisotropic ansatz inherits the same bug class. Wave E Round 2 "
        "discovered the SNStreamingOperator.apply path uses the "
        "reflective-BC equation map, so it cannot close ERR-026 on the "
        "vacuum-BC MMS case without an equation-map extension. Round 3 "
        "owns the closure mechanism."
    ),
)
@pytest.mark.verifies(
    "transport-spherical",
    "sn-mms-spherical-aniso-psi",
    "sn-mms-spherical-aniso-qext",
)
def test_sn_spherical_aniso_mms_converges_second_order():
    r"""Spherical SN with anisotropic ansatz must show :math:`\mathcal{O}(h^2)`.

    Activates the :math:`(1-\mu^2) B(r)/r` angular-redistribution term
    via :math:`\psi_n(r) = (A(r) + B(r) \mu_n)/W` with
    :math:`A(r) = \sin(\pi r/R)`, :math:`B(r) = (r/R)(1-r/R) \cos(\pi r/R)`.
    Both A and B vanish at :math:`r=0` (symmetry) and :math:`r=R`
    (vacuum), so the BCs hold automatically per ordinate.

    The reference scalar flux :math:`\phi(r) = A(r)` is identical to the
    isotropic case (the :math:`\mu_n B(r)` term integrates to zero
    under any quadrature symmetric in :math:`\mu`). The test still
    pins :math:`\phi`, but the underlying angular flux carries the
    redistribution-coupling stress.
    """
    case = build_spherical_anisotropic_mms_case()

    n_cells = [10, 20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux[:, 0, 0]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2), got orders={orders}, errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3


@pytest.mark.xfail(
    strict=True,
    reason=(
        "ERR-026 — same root cause as the spherical anisotropic case; "
        "cylindrical sweep also produces order ≈ 0 on the legacy "
        "isotropic test. Wave E Round 2 discovered the FD operator's "
        "equation map is reflective-BC only — Round 3 owns the closure."
    ),
)
@pytest.mark.verifies(
    "transport-cylindrical",
    "sn-mms-cylindrical-aniso-psi",
    "sn-mms-cylindrical-aniso-qext",
)
def test_sn_cylindrical_aniso_mms_converges_second_order():
    r"""Cylindrical SN with anisotropic ansatz must show :math:`\mathcal{O}(h^2)`.

    Activates the :math:`\xi_n^2 B(r)/r` cylindrical analog of the
    spherical redistribution term. Per the method-implementer's
    derivation, :math:`\xi^2 \neq 1 - \eta^2` in general (since
    :math:`\eta^2 + \xi^2 + \mu^2 = 1`), so the cylindrical case
    exercises a structurally distinct quadrature evaluation from the
    sphere even though the bulk operator structure is shared
    (connection-coefficient analog viewed in two coordinate charts on
    SO(3) — see method-implementer closeout memo).
    """
    case = build_cylindrical_anisotropic_mms_case()

    n_cells = [10, 20, 40, 80]
    errors = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=500, inner_tol=1e-13,
        )
        phi_num = result.scalar_flux[:, 0, 0]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))

    errors = np.asarray(errors)
    orders = np.log2(errors[:-1] / errors[1:])

    assert np.all(orders > 1.9), (
        f"Expected O(h^2), got orders={orders}, errors={errors}"
    )
    assert 1e-8 < errors[-1] < 1e-3
