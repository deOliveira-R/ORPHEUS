"""The Peierls Nystrom cross-checks of the characteristic reference, withdrawn under #506 with their generators.

P1 step (e1b) of the characteristic-reference campaign: the old trajectory-resolvent
family's three #506 rows, re-pointed at the characteristic reference (the roster's
ruling: Peierls-Nystrom stays an independent check, withdrawn under #506) or
re-homed unchanged where they never read the family. Every row here is withdrawn,
so ``tests/gates/test_withdrawal.py`` M5, which runs every file holding a #506 row,
runs this one in a second; the listed ids are in ``tests/gates/withdrawal_506_placement.txt``.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.continuous.characteristic import CharacteristicDerivation, Resolution
from orpheus.geometry import BC, StructuredGeometry
from orpheus.numerics.observable import Eigenvalue
from orpheus.numerics.question import Eigen
from orpheus.numerics.traced_memo import bypass
from orpheus.specification.specification import GeometrySpecification
from tests._harness.withdrawals import PEIERLS_NYSTROM_WITHDRAWN
from tests.gates.derivations._characteristic_ladders import rung
from tests.gates.derivations.test_characteristic_system import _mixture



_K = CellCoefficient.every(Channel.FISSION_EMISSION)


@lru_cache(maxsize=None)
def _derivation(specification: GeometrySpecification, resolution: Resolution) -> CharacteristicDerivation:
    return CharacteristicDerivation(specification, resolution)


def _read(specification: GeometrySpecification, observable, resolution: Resolution) -> float:
    with bypass():
        return float(_derivation(specification, resolution).evaluate(observable).value)


def _one_group(sigma_t: float, sigma_s: float, nu_sigma_f: float):
    return _mixture([sigma_t], [[sigma_s]], None, [nu_sigma_f], [1.0])


@PEIERLS_NYSTROM_WITHDRAWN
@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_the_vacuum_slabs_k_is_the_nystrom_slabs() -> None:
    """[X, withdrawn under #506] The one-group vacuum slab [0, 10] cm (Sigma_t 0.5, Sigma_s 0.38, nu Sigma_f 0.025)
    reads the Peierls Nystrom slab's k (E_1 kernel, mpmath, 8 panels of order 6 at 20 digits) to 5e-5.

    Re-pointed from the old family's ``test_alpha_zero_vacuum_agrees_with_nystrom_slab``
    (the roster's ruling: Peierls-Nystrom stays an independent check, withdrawn
    with its generators under #506).
    """
    from orpheus.derivations.continuous.peierls_nystrom.slab import solve_peierls_eigenvalue

    sigma_t, sigma_s, nu_sigma_f, width = 0.5, 0.38, 0.025, 10.0
    spec = GeometrySpecification(Materials({0: _one_group(sigma_t, sigma_s, nu_sigma_f)}),
                                 StructuredGeometry.slab((0.0, width), (0,), left=BC.vacuum, right=BC.vacuum), Eigen(_K))
    k = _read(spec, Eigenvalue(), rung(3))
    nystrom = solve_peierls_eigenvalue(
        [np.array([sigma_t])], [np.array([[sigma_s]])], [np.array([nu_sigma_f])], [np.array([1.0])], [width],
        n_panels_per_region=8, p_order=6, precision_digits=20, boundary="vacuum",
    )
    k_nystrom = nystrom.k_eff
    assert k_nystrom is not None
    assert abs(k_nystrom - k) / k_nystrom < 5e-5, (k, k_nystrom)


def _thin_sphere():
    return {"sig_t": np.array([0.5]), "sig_s": np.array([[0.38]]), "nu_sig_f": np.array([0.025]), "radii": np.array([5.0])}


@PEIERLS_NYSTROM_WITHDRAWN
@pytest.mark.foundation
def test_the_nystrom_rank_one_specular_closure_is_the_white_hebert_closure() -> None:
    """[withdrawn under #506] Re-homed unchanged from the old family's file (it never read the family):
    the Nystrom ``specular_multibounce`` closure at rank 1 equals ``white_hebert`` to 1e-8, since at N = 1
    R = [[1]] and T_00 = P_ss (V_alpha2)."""
    from orpheus.derivations.continuous.peierls_nystrom.geometry import SPHERE_1D, solve_peierls_1g

    fix = _thin_sphere()
    sol_mb = solve_peierls_1g(SPHERE_1D, fix["radii"], fix["sig_t"], fix["sig_s"], fix["nu_sig_f"],
                              boundary="specular_multibounce", n_bc_modes=1, p_order=4, n_panels_per_region=2,
                              n_angular=24, n_rho=24, n_surf_quad=24, dps=20, tol=1e-10)
    sol_heb = solve_peierls_1g(SPHERE_1D, fix["radii"], fix["sig_t"], fix["sig_s"], fix["nu_sig_f"],
                               boundary="white_hebert", n_bc_modes=1, p_order=4, n_panels_per_region=2,
                               n_angular=24, n_rho=24, n_surf_quad=24, dps=20, tol=1e-10)
    assert sol_mb.k_eff is not None and sol_heb.k_eff is not None
    np.testing.assert_allclose(sol_mb.k_eff, sol_heb.k_eff, rtol=1e-8, atol=1e-10)


@PEIERLS_NYSTROM_WITHDRAWN
@pytest.mark.foundation
def test_the_nystrom_specular_closure_converges_in_its_rank_toward_the_characteristic_closed_sphere() -> None:
    """[withdrawn under #506] The Nystrom ``specular_multibounce`` closure's error against the closed sphere's k
    (the characteristic reference, k = k_inf, D5) is below 0.5 % at rank 1 and 0.2 % at rank 3, and falls
    from rank 1 to 3. Re-pointed from ``test_b5_phase4_converges_toward_variant_alpha``."""
    from orpheus.derivations.continuous.peierls_nystrom.geometry import SPHERE_1D, solve_peierls_1g

    fix = _thin_sphere()
    spec = GeometrySpecification(
        Materials({0: _one_group(0.5, 0.38, 0.025)}),
        StructuredGeometry.sphere((0.0, 5.0), (0,), outer=BC.reflective), Eigen(_K))
    k_closed = _read(spec, Eigenvalue(), rung(3))
    errors = {}
    for rank in (1, 3):
        solution = solve_peierls_1g(SPHERE_1D, fix["radii"], fix["sig_t"], fix["sig_s"], fix["nu_sig_f"],
                                    boundary="specular_multibounce", n_bc_modes=rank, p_order=4, n_panels_per_region=2,
                                    n_angular=24, n_rho=24, n_surf_quad=24, dps=20, tol=1e-10)
        assert solution.k_eff is not None
        errors[rank] = abs(solution.k_eff - k_closed) / k_closed
    assert errors[1] < 5e-3 and errors[3] < 2e-3 and errors[3] <= errors[1], errors
