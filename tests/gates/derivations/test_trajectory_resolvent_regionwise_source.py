r"""The multi-region trajectory-resolvent oracles read the emission density region by region (ERR-090).

In a piecewise-homogeneous medium the isotropic emission density
:math:`q(r) = \Sigma_s(r)\phi(r) + \chi(r)\nu\Sigma_f(r)\phi(r)/k` is smooth
inside each region and jumps at every material interface. The multi-region
sphere and cylinder chord oracles
(:class:`~orpheus.derivations.continuous.trajectory_resolvent.chord_oracle.MultiRegionSphereChordOracle`,
:class:`~orpheus.derivations.continuous.trajectory_resolvent.chord_oracle.MultiRegionCylinderChordOracle`)
reconstruct :math:`q` along a chord from its values at the radial nodes. Until
ERR-090 they fitted ONE cubic spline through every node, across the jumps, and
the references converged in :math:`n_r` at about first order and
non-monotonically. They now build one spline per region
(:func:`~orpheus.derivations.continuous.trajectory_resolvent.chord_oracle._regionwise_cubic_spline`),
and each chord segment evaluates its own region's piece.

The ladder in this file, bottom up:

1. ``test_regionwise_interpolant_reproduces_a_piecewise_cubic`` (THEOREM): a
   not-a-knot cubic spline through four or more nodes reproduces a cubic
   exactly, so the per-region interpolant of a density that is a different
   cubic in each region is exact at every radius, the interface
   neighbourhoods included. Its refusal leg pins the two-node floor.
2. ``test_mr_oracle_first_leg_matches_the_line_integral`` (REFERENCE, one row
   per geometry): with the reflectivity at zero the oracle's output is the
   first-leg line integral of :math:`q` against the attenuation. For a
   piecewise-cubic :math:`q` that integral is computed here independently,
   with :func:`scipy.integrate.quad` along a chord parametrised from first
   principles, split at the interface crossings. The only thing the two
   sides share is the density's definition.

The :math:`n_r` ladder of the sphere's eigenvalue, the observed convergence
this defect broke, is ``test_mr_sphere_k_converges_in_n_r`` in
``test_peierls_greens_function_mr.py``; it rests on row 2.

Mutation witness (``[M]`` 2026-09-26, ``python -O -m pytest`` with the global
spline re-installed in process, every piece replaced by one spline through all
nodes): rows 1 and 2 red at relative errors of 5.3e-2 (sphere) and 6.7e-2
(cylinder) against an honest 3.0e-10 / 8.0e-10.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import quad

from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    MultiRegionCylinderChordOracle,
    MultiRegionSphereChordOracle,
    _regionwise_cubic_spline,
)
from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
    _composite_per_region_gl,
)

_THIS = "tests/gates/derivations/test_trajectory_resolvent_regionwise_source.py"

# Three regions with a jump at each interface; the density is a different cubic in each.
_RADII = np.array([0.5, 1.5, 2.0])
_R = float(_RADII[-1])
_CUBIC_PER_REGION = np.array([
    [1.00, 0.30, -0.40, 0.20],
    [0.20, -0.10, 0.25, -0.05],
    [1.50, 0.40, -0.30, 0.10],
])  # rows: region; columns: coefficients of 1, r, r^2, r^3
_SIGMA_T = 0.7  # [1/cm], uniform: the attenuation is not what is under test

#: The oracle integrates each segment with a 64-point Gauss-Legendre rule. On a
#: chord passing close to the origin the density, a function of
#: r(s) = sqrt(h^2 + (s - s0)^2), is nearly kinked, and that rule's error there
#: is the floor: [M] 2026-09-26, 3.0e-10 (sphere) and 8.0e-10 (cylinder),
#: relative to the largest angular flux. The tolerance is ten times the larger.
_LINE_INTEGRAL_RTOL = 1e-8


def _density(r: float) -> float:
    """The piecewise-cubic emission density, region k owning (r_{k-1}, r_k]."""
    k = min(int(np.searchsorted(_RADII, r, side="left")), len(_RADII) - 1)
    return float(np.polynomial.polynomial.polyval(r, _CUBIC_PER_REGION[k]))


def _nodes() -> tuple[np.ndarray, np.ndarray]:
    """The composite per-region Gauss-Legendre nodes the solvers use (6/12/6 here), with their regions."""
    r_nodes, _, region_at_node = _composite_per_region_gl(_RADII, 24)
    return r_nodes, region_at_node


def _first_leg_line_integral(r: float, cos_theta: float, inv_sin_axial: float) -> float:
    r"""The first-leg integral, from first principles.

    A ray at radius ``r`` travelling (in the plane of the chord) at angle
    :math:`\theta` to the outward radial direction came from the point
    :math:`s` behind it, at radius
    :math:`\rho(s) = \sqrt{r^2 - 2 r s\cos\theta + s^2}`, and entered at
    :math:`L = r\cos\theta + \sqrt{R^2 - r^2\sin^2\theta}`. With in-plane
    arclength ``s`` and the axial lift ``inv_sin_axial`` (1 for a sphere),

    .. math::

       \psi = \text{inv\_sin\_axial} \int_0^L q(\rho(s))\,
              e^{-\Sigma_t\, s\, \text{inv\_sin\_axial}}\, ds,

    split at every interface crossing :math:`\rho(s) = r_k`.
    """
    L = r * cos_theta + np.sqrt(_R * _R - r * r * (1.0 - cos_theta * cos_theta))
    crossings = []
    for r_k in _RADII[:-1]:
        disc = (r * cos_theta) ** 2 - (r * r - r_k * r_k)
        if disc > 0.0:
            crossings += [r * cos_theta - np.sqrt(disc), r * cos_theta + np.sqrt(disc)]
    inside = [s for s in crossings if 0.0 < s < L]

    def integrand(s: float) -> float:
        rho = np.sqrt(max(r * r - 2.0 * r * cos_theta * s + s * s, 0.0))
        return _density(rho) * np.exp(-_SIGMA_T * s * inv_sin_axial)

    value, _ = quad(integrand, 0.0, L, points=inside or None, epsabs=1e-14, epsrel=1e-13, limit=200)
    return inv_sin_axial * value


@pytest.mark.foundation
def test_regionwise_interpolant_reproduces_a_piecewise_cubic() -> None:
    r"""One not-a-knot spline per region reproduces a per-region cubic exactly.

    Evaluated off the nodes, at 200 radii per region reaching each region's
    end points (the extrapolated gap between the outermost node and the
    interface included). A single spline across the jumps misses by the size
    of the jump there. The refusal leg: a region holding one node has no
    spline, and the message names the region.
    """
    r_nodes, region_at_node = _nodes()
    values = np.array([_density(r) for r in r_nodes])
    pieces = _regionwise_cubic_spline(r_nodes, values, region_at_node, len(_RADII))
    assert len(pieces) == len(_RADII)
    inner = np.concatenate([[0.0], _RADII[:-1]])
    worst = 0.0
    for k, (a, b) in enumerate(zip(inner, _RADII)):
        r = np.linspace(a, b, 200)
        exact = np.polynomial.polynomial.polyval(r, _CUBIC_PER_REGION[k])
        worst = max(worst, float(np.max(np.abs(pieces[k](r) - exact))))
    # A cubic reproduced from values at 6 to 12 nodes: rounding only.
    assert worst < 1e-12, f"the regionwise interpolant misses a per-region cubic by {worst:.3e}"

    with pytest.raises(ValueError, match=r"region 0 holds 1 radial node"):
        _regionwise_cubic_spline(
            np.array([0.25, 1.0, 1.2, 1.7, 1.8]), np.ones(5), np.array([0, 1, 1, 2, 2]), len(_RADII),
        )


@pytest.mark.l0
@pytest.mark.catches("ERR-090")
@pytest.mark.verifies("peierls-greens-mr-regionwise-source")
@pytest.mark.rests_on(f"{_THIS}::test_regionwise_interpolant_reproduces_a_piecewise_cubic")
@pytest.mark.parametrize("geometry", ["sphere", "cylinder"])
def test_mr_oracle_first_leg_matches_the_line_integral(geometry: str) -> None:
    r"""The oracle at zero reflectivity is the first-leg line integral of the density.

    Activates: the per-segment choice of the density's piece (the density
    jumps at both interfaces, and chords through the moderator cross both),
    the segment split at the interfaces, and, for the cylinder, the axial lift.
    Nulls: the bounce period (zero reflectivity) and any jump in the total
    cross section (uniform here), both exercised by the solver-level rows.

    The tolerance is the oracle's own chord-quadrature floor times ten
    (:data:`_LINE_INTEGRAL_RTOL`); with the global spline the reading is
    5.3e-2 (sphere) and 6.7e-2 (cylinder).
    """
    r_nodes, region_at_node = _nodes()
    density_at_nodes = np.array([_density(r) for r in r_nodes])
    sigma_t = np.full(len(_RADII), _SIGMA_T)
    if geometry == "sphere":
        mu = np.array([-0.9, -0.4, 0.1, 0.6, 0.95])
        oracle = MultiRegionSphereChordOracle(
            r_nodes=r_nodes, mu_nodes=mu, R=_R, radii=_RADII,
            sigma_t_per_region=sigma_t, alpha=0.0, region_at_node=region_at_node,
        )
        reference = np.array([[_first_leg_line_integral(r, m, 1.0) for m in mu] for r in r_nodes])
    else:
        mu_axial = np.array([-0.7, 0.2, 0.9])
        phi_az = np.array([0.3, 1.2, 2.0, 3.5, 5.9])
        oracle = MultiRegionCylinderChordOracle(
            r_nodes=r_nodes, mu_axial_nodes=mu_axial, phi_az_nodes=phi_az, R=_R,
            radii=_RADII, sigma_t_per_region=sigma_t, alpha=0.0, region_at_node=region_at_node,
        )
        reference = np.array([
            [[_first_leg_line_integral(r, np.cos(p), 1.0 / np.sqrt(1.0 - m * m)) for p in phi_az]
             for m in mu_axial]
            for r in r_nodes
        ])
    produced = oracle.apply_operator(density_at_nodes, sigma_t=0.0, n_traj_quad=64)
    error = float(np.max(np.abs(produced - reference)) / np.max(np.abs(reference)))
    print(f"{geometry}: max |oracle - line integral| / max |line integral| = {error:.3e}")
    assert error < _LINE_INTEGRAL_RTOL, (
        f"the {geometry} multi-region oracle's first leg misses the line integral of a "
        f"piecewise-cubic density by {error:.3e} (relative; tolerance {_LINE_INTEGRAL_RTOL:.0e}): "
        f"the emission density is not being read from each segment's own region"
    )
