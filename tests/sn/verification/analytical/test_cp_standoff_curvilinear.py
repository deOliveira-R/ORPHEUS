"""Curvilinear SN ↔ CP cross-checks (ancillary code-to-code standoff).

Split from the legacy ``tests/sn/test_cylindrical.py`` /
``test_spherical.py`` (SN taxonomy reorg). These are the SN-vs-CP
cross-checks: homogeneous-1G (both must hit the analytical k_inf, so
the white-BC vs reflective-BC gap vanishes) and heterogeneous (~10%
agreement — the white-BC approximation gap). Per ``vv-principles``
this is L4-style code-to-code agreement *backed by* the analytical
k_inf at the homogeneous limit; the heterogeneous leg is a benchmark
band, not a verification target.

The cylindrical legs are below; the spherical leg is appended by the
test_spherical split (phase 4).
"""

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn
from tests.sn._test_helpers import (
    curvilinear_homogeneous_mesh as _homogeneous_mesh,
    curvilinear_two_region_mesh as _two_region_mesh,
)

# Equation-coverage list preserved verbatim from the legacy
# legacy modules so no verifies(...) edge is lost in the split. Cylinder
# and sphere carry DIFFERENT lists, so they are applied per-section.
_CYL_VERIFIES = pytest.mark.verifies(
    "transport-cylindrical",
    "alpha-cylindrical",
    "alpha-recursion",
    "wdd-closure",
    "wdd-face",
    "mm-weights",
    "multigroup",
    "reflective-bc",
    "one-group-kinf",
    "matrix-eigenvalue",
    "mg-balance",
    "balance-general",
)
_SPH_VERIFIES = pytest.mark.verifies(
    "transport-spherical",
    "alpha-recursion",
    "wdd-closure",
    "wdd-face",
    "multigroup",
    "reflective-bc",
    "one-group-kinf",
    "matrix-eigenvalue",
    "mg-balance",
    "balance-general",
)


# ═══════════════════════════════════════════════════════════════════════
# Cylindrical SN ↔ CP
# ═══════════════════════════════════════════════════════════════════════

@_CYL_VERIFIES
@pytest.mark.l1
def test_cross_check_with_cp_1g():
    """SN and CP on the same cylindrical geometry should give close k_inf."""
    from orpheus.cp.solver import solve_cp

    mix = get_mixture("A", "1g")

    mesh_sn = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL)
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    result_sn = solve_sn({0: mix}, mesh_sn, quad,
                         max_inner=500, inner_tol=1e-10)

    mesh_cp = _homogeneous_mesh(
        1, 2.0, mat_id=0, coord=CoordSystem.CYLINDRICAL, bc=BC.white,
    )
    result_cp = solve_cp({0: mix}, mesh_cp)

    np.testing.assert_allclose(
        result_sn.keff, result_cp.keff, rtol=1e-6,
        err_msg=f"SN keff={result_sn.keff:.6f} vs CP keff={result_cp.keff:.6f}",
    )


@_CYL_VERIFIES
@pytest.mark.l2
@pytest.mark.slow
def test_heterogeneous_sn_vs_cp_cross_check():
    """Heterogeneous SN and CP should agree within ~10%."""
    from orpheus.cp.solver import solve_cp

    mix_fuel = get_mixture("A", "1g")
    mix_mod = get_mixture("B", "1g")
    materials = {2: mix_fuel, 0: mix_mod}

    mesh_sn = _two_region_mesh(
        outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(20, 20),
        coord=CoordSystem.CYLINDRICAL,
    )
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    result_sn = solve_sn(materials, mesh_sn, quad,
                         max_inner=500, inner_tol=1e-10)

    mesh_cp = _two_region_mesh(
        outers=(0.5, 1.0), mat_ids=(2, 0), n_cells=(10, 10),
        coord=CoordSystem.CYLINDRICAL,
        bc=BC.white,  # CP supports vacuum/white only
    )
    result_cp = solve_cp(materials, mesh_cp)

    np.testing.assert_allclose(
        result_sn.keff, result_cp.keff, rtol=0.10,
        err_msg=f"SN={result_sn.keff:.6f} vs CP={result_cp.keff:.6f}",
    )


# ═══════════════════════════════════════════════════════════════════════
# Spherical SN ↔ CP (split from the legacy test_spherical.py)
# ═══════════════════════════════════════════════════════════════════════

@_SPH_VERIFIES
@pytest.mark.l1
def test_cross_check_with_cp_1g_sphere():
    """SN and CP on the same spherical geometry should give close k_inf.

    For homogeneous 1G, both must match the analytical value exactly.
    The white-BC (CP) vs reflective-BC (SN) difference vanishes for
    homogeneous infinite medium.
    """
    from orpheus.cp.solver import solve_cp

    mix = get_mixture("A", "1g")

    mesh_sn = _homogeneous_mesh(20, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
    quad = Quadrature.gauss_legendre(8)
    result_sn = solve_sn({0: mix}, mesh_sn, quad,
                         max_inner=500, inner_tol=1e-10)

    mesh_cp = _homogeneous_mesh(
        1, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL, bc=BC.white,
    )
    result_cp = solve_cp({0: mix}, mesh_cp)

    np.testing.assert_allclose(
        result_sn.keff, result_cp.keff, rtol=1e-6,
        err_msg=f"SN keff={result_sn.keff:.6f} vs CP keff={result_cp.keff:.6f}",
    )
