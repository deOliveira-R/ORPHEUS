r"""L1 MMS convergence tests for curvilinear SN with **anisotropic** ansatz.

Companion to :file:`tests/sn/test_mms_curvilinear.py` (legacy isotropic
ansatz). This file exercises the angular-redistribution operator
directly via the :math:`(A(r) + B(r) \cdot \zeta) / W` ansatz from the
method-implementer's Phase-1 work, where :math:`\zeta = \mu` for sphere
and :math:`\zeta = \eta` for cylinder.

Failure-mode #7 from :ref:`vv-principles` motivates the redesign: the
isotropic ansatz nulls the angular redistribution by construction, so
the MMS test cannot detect ERR-026-class bugs (curvilinear sweep WDD
wrong fixed point). The anisotropic case adds the :math:`B(r) \zeta`
term that activates :math:`(1-\mu^2) B/r` (sphere) or :math:`\xi^2 B/r`
(cylinder) — the angular-redistribution analytic — making the
curvilinear sweep's hardest math observable under refinement.

Wave E Round 3 status (the ERR-026 closure attempt): Round 3 shipped
the BC-aware FD operator (:func:`solution_to_angular_flux*` now consume
``sn_mesh.bc_*`` via :meth:`ResolvedBC.apply_to_incoming`), which is
the infrastructure load-bearing for any closure of ERR-026 on
fixed-source MMS — vacuum, reflective, white, albedo, and mixed BCs
are now plumbed uniformly through the FD operator. **However**, the
curvilinear-default flip from ``"source_iteration"`` to ``"krylov"``
in :func:`solve_sn_fixed_source` was not enabled because empirically
the symmetric-closure FD operator at the curvilinear outer face uses
cell-center as a face-flux approximation, which is only first-order
accurate on non-constant solutions.  Switching the default to
``"krylov"`` regresses the MMS convergence rate from the WDD sweep's
~:math:`\mathcal{O}(h^{1.3})` (ERR-026-affected, but a benign
volumetric error mode for these MMS) to ~:math:`\mathcal{O}(h^{1})`
(FD operator's boundary truncation).

Both tests therefore stay ``@pytest.mark.xfail`` — Round 3's BC
plumbing is the infrastructure that *enables* a future full closure
of ERR-026 on MMS; the closure itself depends on a follow-up that
fixes the FD operator's outer-face flux treatment (DD diamond
relation at the boundary, or extrapolation that preserves second-
order accuracy at the boundary).

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
        "ERR-026 — curvilinear sweep WDD wrong fixed point + curvilinear "
        "FD operator boundary face-flux approximation is first-order "
        "accurate on non-constant solutions.  Wave E Round 3 shipped "
        "the BC-aware FD operator (load-bearing infrastructure: "
        "solution_to_angular_flux* + the matvec helpers consume "
        "sn_mesh.bc_* via ResolvedBC.apply_to_incoming), which fixes "
        "the vacuum-BC slot fill that Round 2 identified as broken; "
        "but the FD operator's cell-center-as-face-value approximation "
        "at the outer boundary regresses the MMS convergence rate "
        "from the sweep's ~O(h^1.3) (ERR-026 benign mode) to ~O(h^1) "
        "if used as the default.  Round 3 therefore keeps "
        "``inner_solver=\"source_iteration\"`` as the default for all "
        "geometries; ``\"krylov\"`` is opt-in and correct for "
        "constant-source problems but not the right default for MMS. "
        "Full ERR-026 closure on MMS depends on a follow-up that "
        "extrapolates the curvilinear outer-face flux at second order."
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
        "ERR-026 — same root cause as the spherical anisotropic case. "
        "Wave E Round 3 shipped the BC-aware FD operator infrastructure; "
        "full ERR-026 MMS closure depends on a follow-up that "
        "extrapolates the curvilinear outer-face flux at second order."
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
