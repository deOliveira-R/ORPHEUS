r"""``compute_P_ss_sphere`` against its closed form: the sphere sibling of
``tests/gates/cp/test_cylinder_pss.py``.

:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.compute_P_ss_sphere`
is the surface-to-surface uncollided transmission of a sphere under a
uniform isotropic inward partial current, the building block of the
sphere ``white_hebert`` closure (Hébert 2009 §3.8.5 Eq. (3.323)).  For a
homogeneous sphere, with :math:`\tau_R = \Sigma_t R`,

.. math::

   P_{ss} = 2\int_0^1 \mu\,e^{-2\tau_R\mu}\,d\mu
          = \frac{1 - (1 + 2\tau_R)\,e^{-2\tau_R}}{2\,\tau_R^{2}}

(:eq:`peierls-class-b-Pss-homogeneous` on the theory page;
``peierls-sphere-Pss-homogeneous`` in the function's docstring).  The
integral-to-closed-form step is pinned symbolically by
``tests/gates/derivations/test_peierls_greens_function_symbolic.py``
(V_alpha2 Path B); this module pins PRODUCTION to it.

Before this module the only gate on ``compute_P_ss_sphere`` compared it
with ``compute_T_specular_sphere``, which evaluates the same impact-parameter
integral over the same chord geometry (``chord_quadrature`` and
``chord_half_lengths``), so the two agree by construction (X4): an error in
the shared chord geometry cancels there.

Structural independence of the reference: the closed form is evaluated in
:mod:`mpmath` at 40 digits from :math:`\tau_R` alone; production integrates
:math:`(2/R^2)\int_0^R h\,e^{-\tau(h)}\,dh` by Gauss-Legendre in the impact
parameter :math:`h`, with :math:`\tau(h)` summed from per-annulus chord
segments.

Origin: ``derivations/diagnostics/diag_sphere_white_bc_geometric_series_fix.py``
(retired, R19; #132), whose ``P_ss`` closed-form leg is promoted here (the
closure it prototyped shipped as ``white_hebert`` and is gated by the
k_eff tests); the probe could not be ported (its imports were dead).
"""

from __future__ import annotations

import mpmath
import numpy as np
import pytest

from orpheus.derivations.continuous.peierls_nystrom.geometry import (
    compute_P_ss_sphere,
)

pytestmark = [
    pytest.mark.l0,
    pytest.mark.verifies(
        "peierls-class-b-Pss-homogeneous", "peierls-sphere-Pss-homogeneous",
    ),
]

_EPS = float(np.finfo(float).eps)


def _pss_closed_form(tau_R: float) -> float:
    with mpmath.workdps(40):
        t = mpmath.mpf(tau_R)
        return float((1 - (1 + 2 * t) * mpmath.exp(-2 * t)) / (2 * t * t))


def _pss_tolerance(n_quad: int, tau_R: float) -> float:
    r"""``10 eps (n_quad + 4 tau_R (1 + tau_R))``.

    Two rounding sources.  The Gauss-Legendre sum has depth ``n_quad``.
    The chord :math:`2\sqrt{R^2 - h^2}` carries a relative rounding of
    :math:`\epsilon R^2 / (R^2 - h^2)`, largest at the rim, and
    :math:`e^{-\tau}` turns it into a relative error :math:`\tau` times
    that; averaged over the integrand :math:`h\,e^{-\tau(h)}` the
    amplification is :math:`\kappa = 4\tau_R` for a thin cell and
    :math:`4\tau_R^2` for a thick one, both below
    :math:`4\tau_R(1 + \tau_R)` (`[M]` 2026-09-22 by :func:`mpmath.quad`
    of the weighted mean: 0.41, 5.8, 1600 and 40000 at
    :math:`\tau_R` = 0.1, 1, 20, 100).  A thick cell's transmission is
    carried by rim chords, which is why it is the ill-conditioned end.
    """
    return 10.0 * _EPS * (n_quad + 4.0 * tau_R * (1.0 + tau_R))


_SPHERE_PSS_RULES = [
    "tests/gates/derivations/test_peierls_greens_function_symbolic.py::test_v_alpha2_Pss_matches_hebert_via_polar_path",
    "tests/gates/derivations/test_quadrature.py::test_chord_quadrature_recovers_constant_integral",
    "tests/gates/derivations/test_quadrature.py::test_chord_quadrature_homogeneous_panel_structure",
]


@pytest.mark.rests_on(*_SPHERE_PSS_RULES)
@pytest.mark.parametrize("R", [1.0, 2.5])
@pytest.mark.parametrize("tau_R", [1e-4, 1e-2, 0.1, 0.5, 1.0, 2.0, 5.0, 20.0])
def test_homogeneous_pss_matches_the_closed_form(tau_R: float, R: float) -> None:
    r"""Production equals the closed form from the thin (:math:`\tau_R =
    10^{-4}`, :math:`P_{ss} \to 1`) to the thick (:math:`\tau_R = 20`)
    cell, at ``n_quad = 64``, within :func:`_pss_tolerance`.

    What it catches: an error in the chord length (a half chord where the
    full antipodal chord belongs), in the :math:`h` weight of the
    impact-parameter integrand, or in the :math:`2/R^2` normalisation.  The
    closed form depends on :math:`\tau_R` only, so the :math:`R = 2.5`
    rows (with :math:`\Sigma_t = \tau_R / R`) are the ones that see an
    :math:`R`-power error; at :math:`R = 1` every power of :math:`R` is 1.
    `[M]` 2026-09-22: worst relative error :math:`7.4\times10^{-14}`
    over the 16 rows.
    """
    n = 64
    p = compute_P_ss_sphere(np.array([R]), np.array([tau_R / R]), n_quad=n)
    ref = _pss_closed_form(tau_R)
    rel = abs(p - ref) / ref
    tol = _pss_tolerance(n, tau_R)
    assert rel < tol, (
        f"P_ss(tau_R={tau_R}, R={R}, n_quad={n}) = {p:.16e}, closed form "
        f"{ref:.16e}: rel {rel:.3e} >= {tol:.3e}"
    )


@pytest.mark.rests_on(
    "tests/gates/cp/test_sphere_pss.py::test_homogeneous_pss_matches_the_closed_form",
)
@pytest.mark.parametrize("R", [1.0, 2.5])
def test_an_optically_thick_cell_matches_once_resolved(R: float) -> None:
    r""":math:`\tau_R = 100`: the transmission is carried by a boundary layer
    of width :math:`1/(2\tau_R)` in :math:`\mu`, so ``n_quad = 64`` does not
    resolve it and 128 does.  `[M]` 2026-09-22 at :math:`R = 1`: relative
    error :math:`7.3\times10^{-13}` at 64, :math:`3.3\times10^{-13}` at 128
    and :math:`4.8\times10^{-13}` at 256 (the rounding floor, not the
    quadrature); :math:`8.4\times10^{-13}` at :math:`R = 2.5`, 128; the
    bound is :math:`9.0\times10^{-11}`, set by the rim conditioning
    :math:`4\tau_R^2` (:func:`_pss_tolerance`).
    """
    tau_R, n = 100.0, 128
    p = compute_P_ss_sphere(np.array([R]), np.array([tau_R / R]), n_quad=n)
    ref = _pss_closed_form(tau_R)
    rel = abs(p - ref) / ref
    tol = _pss_tolerance(n, tau_R)
    assert rel < tol, (
        f"P_ss(tau_R={tau_R}, R={R}, n_quad={n}): rel {rel:.3e} >= {tol:.3e}"
    )


@pytest.mark.rests_on(
    "tests/gates/cp/test_sphere_pss.py::test_homogeneous_pss_matches_the_closed_form",
    "tests/gates/derivations/test_quadrature.py::test_chord_quadrature_multi_region_panel_structure",
)
@pytest.mark.parametrize("tau_R", [0.5, 1.0, 2.0, 20.0])
def test_a_multiregion_sphere_of_one_material_is_the_homogeneous_sphere(
    tau_R: float,
) -> None:
    r"""Splitting a homogeneous sphere into 2 or 4 shells of the same
    :math:`\Sigma_t` leaves :math:`P_{ss}` unchanged: the multi-region
    routing (per-annulus chord segments, the quadrature subdivided at every
    interior radius) is invisible when the material is.

    What it catches: a routing error, a chord segment dropped, double
    counted or attributed to the wrong shell, that the one-region closed
    form cannot see.  Declared blind: an error common to both sides (the
    chord length, the weight, the normalisation), which the closed-form
    test owns.  Tolerance: each side is within :func:`_pss_tolerance` of
    the truth, so they are within twice it of each other.  `[M]` 2026-09-22:
    worst difference :math:`2.4\times10^{-14}` (4 shells, :math:`\tau_R =
    20`).
    """
    n, R, sig = 64, 1.0, np.array([tau_R])
    p1 = compute_P_ss_sphere(np.array([R]), sig, n_quad=n)
    tol = 2.0 * _pss_tolerance(n, tau_R)
    for radii in (np.array([0.5, 1.0]), np.array([0.4, 0.45, 0.55, 1.0])):
        p_mr = compute_P_ss_sphere(radii, np.full(len(radii), tau_R), n_quad=n)
        rel = abs(p_mr - p1) / p1
        assert rel < tol, (
            f"tau_R={tau_R}: {len(radii)} shells of one material give "
            f"P_ss={p_mr:.16e} vs one region {p1:.16e}, rel {rel:.3e} >= {tol:.3e}"
        )
