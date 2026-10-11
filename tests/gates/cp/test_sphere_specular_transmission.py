r"""``compute_T_specular_sphere`` against its integral form, every mode pair (#578, P2 group D).

:math:`T_{mn} = 2\int_0^1 \mu\,\tilde P_m(\mu)\,\tilde P_n(\mu)\,e^{-2\tau_R\mu}\,d\mu` for a homogeneous sphere
(``peierls-specular-T-matrix`` on the Peierls page; :math:`\tilde P_n(\mu) = P_n(2\mu - 1)`), evaluated here by
:func:`mpmath.quad` in :math:`\mu`, sharing no impact-parameter rule and no chord with production, which integrates
in the impact parameter :math:`h` with :math:`\mu = \sqrt{1 - h^2/R^2}` (since #578 the half-chord at :math:`R`
over :math:`R`, ``half_chord(R, h, 0) / R``).

Before this file the only pin of the function compared :math:`T_{00}` with :math:`P_{ss}` (both production, one
chord geometry, X4), blind to :math:`\mu` at every :math:`m, n \ge 1`, where :math:`\mu` is all the modes read.
"""
from __future__ import annotations

from typing import Any

import mpmath
import numpy as np
import pytest

from orpheus.derivations.continuous.peierls_nystrom.geometry import compute_T_specular_sphere

_EPS = float(np.finfo(float).eps)


def _tolerance(n_quad: int, tau_R: float) -> float:
    """The P_ss rows' band (``tests/gates/cp/test_sphere_pss.py::_pss_tolerance``): ``10 eps (n_quad + 4 tau_R (1 +
    tau_R))``, the Gauss sum's depth and the rim's amplification by the optical depth."""
    return 10.0 * _EPS * (n_quad + 4.0 * tau_R * (1.0 + tau_R))


@pytest.mark.l0
@pytest.mark.verifies("peierls-specular-T-matrix")
@pytest.mark.rests_on("tests/gates/cp/test_sphere_pss.py::test_homogeneous_pss_matches_the_closed_form")
@pytest.mark.parametrize(("tau_R", "R"), [(0.5, 1.0), (2.0, 2.5), (20.0, 1.0)])
def test_every_mode_pair_of_the_sphere_transmission_is_its_integral(tau_R: float, R: float) -> None:
    """``T_mn``, m, n < 3, within ``_tolerance`` relative of the 30-digit integral (`[M]` 2026-10-10 at the default
    ``n_quad = 64``: 8.7e-16, 2.3e-14, 6.7e-15). First reds: the exit cosine taken as the sine ``h/R``; the cosine
    formed as ``sqrt(1 - h/R)`` (the square dropped).
    """
    T = compute_T_specular_sphere(np.array([R]), np.array([tau_R / R]), 3)
    with mpmath.workdps(30):
        def p(n, u):
            return mpmath.legendre(n, 2 * u - 1)
        def integral(m: int, n: int) -> Any:
            return 2 * mpmath.quad(lambda u: u * p(m, u) * p(n, u) * mpmath.exp(-2 * tau_R * u), [0, 1])
        ref = np.array([[float(integral(m, n)) for n in range(3)] for m in range(3)])

    rel = np.abs(T - ref) / np.abs(ref)
    assert np.all(rel <= _tolerance(64, tau_R)), f"tau_R = {tau_R}: relative errors\n{rel}"
