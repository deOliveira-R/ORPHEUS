r"""Plan-(b) Option-1 — L1 flux-shape cross-check vs Garcia 2021
J. Comp. Phys. 433.

Garcia 2021's stable-P_N solver provides the only published high-
precision multi-region sphere reference. The paper covers
**fixed-source only** — k-eigenvalue is explicitly out of scope (the
paper's §III calls criticality "future work"). For Variant α multi-
region verification, this gives **flux-shape L1 evidence**: agreement
on the converged scalar flux profile :math:`\phi(r)` across a 3-region
sphere with non-trivial XS jumps, vacuum BC, and constant-per-region
isotropic external source.

**Convention conversion** (memo
``ps1982_and_garcia_extraction.md`` Garcia §):

- Garcia's source ``S_k`` is **per cm³ per steradian** (standard).
  My ``external_source`` argument is **total per cm³** (gets divided
  by :math:`4\pi` internally for isotropic per-steradian source).
  These differ by :math:`4\pi`.
- Garcia's "scalar flux" :math:`\phi(r) = \int_{-1}^{1} \Psi(r,\mu)\,
  \mathrm d\mu` (no :math:`2\pi`). My output
  :math:`\phi(r) = 2\pi \int_{-1}^{1} \psi\,\mathrm d\mu` (standard
  scalar flux). These differ by :math:`2\pi`.

Net for matching: passing ``external_source = (0.5, 1.0, 1.5)`` and
comparing my output to Garcia's table 5 values:

.. math::

    \phi_{\rm mine} \;=\; \tfrac{2\pi}{4\pi}\,\phi_{\rm Garcia}
                    \;=\; \tfrac{1}{2}\,\phi_{\rm Garcia}.

So the agreement criterion is
:math:`\phi_{\rm mine} \approx 0.5 \cdot \phi_{\rm Garcia}^{\rm table}`.
This factor-of-2 conversion is **the convention map**, not an error.

Garcia 2021 Case 1 setup (Williams 1991 Example 5):

- ``radii = [3.0, 5.0, 7.0]`` cm
- Region 1 (core, 0–3 cm): :math:`\Sigma_t = 1.0,\,\Sigma_s = 0.99\,(c=0.99)`
- Region 2 (mid,  3–5 cm): :math:`\Sigma_t = 0.5,\,\Sigma_s = 0.30\,(c=0.6)`
- Region 3 (outer, 5–7 cm): :math:`\Sigma_t = 2.0,\,\Sigma_s = 1.90\,(c=0.95)`
- Internal sources :math:`S = (0.5, 1.0, 1.5)` per cm³ per steradian
- Vacuum BC at r=7

Garcia 2021 verified this Case 1 against Williams 1991 (integral-eq
MOC) and Picca-Furfaro-Ganapol 2012 (S_N) to 3-4 sig figs at every
r-point — three structurally-independent methods agreeing.

Tolerances, derived 2026-09-26 (ERR-090) by
:func:`tests.gates.derivations._trajectory_resolvent_ladders.tolerance_for`
from two measured quantities, neither of them this comparison's reading:

* Garcia's own error: the table carries five significant figures, so each
  value is within half a unit in its last digit; the largest relative
  rounding is 4.6e-5 (on 10.807), computed below from the table's strings.
* The trajectory resolvent's own error at this fixture
  (:math:`n_r = 48`, :math:`n_\mu = 24`), from its OWN ladder: the largest
  relative change to its finer rungs (48, 48), (96, 48), (96, 96) is 1.43e-3
  at the interior table radii (r = 3.5 cm, where the :math:`\mu` rule meets
  the tangency kinks of the two interfaces, #516) and 1.65e-2 at the
  surface, where the flux's radial derivative is singular at the vacuum
  boundary and the outermost region's spline extrapolates to R. Re-measure
  with ``python -O -m tests.gates.derivations._trajectory_resolvent_ladders
  garcia``.

The rule gives 3e-3 inside and 4e-2 at the surface. ``[M]`` the readings are
at most 1.05e-3 inside and 2.3e-2 at the surface; with the one-spline
emission density (ERR-090) they read 1.1e-2 at r = 5.0 cm and 3.3e-3 at
r = 5.5 cm, the second only 10 % over its tolerance.

Until 2026-09-26 the bands were 2 % inside, 15 % within 2 cm of an
interface and 5 % at the surface; the 15 % band covered the one-spline
emission density. The scalar flux is interpolated to the table's radii
per region (its radial derivative jumps at each interface); at an interface
radius the two one-sided values are averaged.
"""
from __future__ import annotations

import numpy as np
import pytest
from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    _regionwise_cubic_spline,
)
from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
    solve_greens_function_sphere_mr_fixed_source,
)
from tests.gates.derivations._trajectory_resolvent_ladders import (
    GARCIA_CASE1_RESOLVENT_STEP,
    tolerance_for,
)


# Garcia 2021 Table 5 (Case 1 converged ppP_N; rightmost column).
GARCIA_2021_CASE1_R = np.array([
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0,
    3.5, 4.0, 4.5, 5.0,
    5.5, 6.0, 6.5, 7.0,
])
#: The table's values as printed (the strings carry the precision).
GARCIA_2021_CASE1_PHI_PRINTED = (
    "18.860", "18.756", "18.442", "17.911", "17.145", "16.095", "14.381",  # core
    "13.455", "13.337", "13.590", "14.361",                                # mid
    "15.532", "14.198", "10.807", "4.0763",                                # outer
)
GARCIA_2021_CASE1_PHI = np.array([float(v) for v in GARCIA_2021_CASE1_PHI_PRINTED])
#: Half a unit in each value's last printed digit, relative: Garcia's own error.
GARCIA_2021_CASE1_ROUNDING = np.array(
    [0.5 * 10.0 ** -len(v.split(".")[1]) for v in GARCIA_2021_CASE1_PHI_PRINTED]
) / GARCIA_2021_CASE1_PHI
# Convention conversion: Variant α output = 0.5 × Garcia table.
GARCIA_TO_VARIANT_ALPHA_FACTOR = 0.5

GARCIA_2021_CASE1_RADII = np.array([3.0, 5.0, 7.0])

#: The two tolerances (module docstring): Garcia's largest rounding as the
#: reference bound, the resolvent's own ladder step as the SUT residual.
INTERIOR_TOLERANCE = tolerance_for(
    GARCIA_CASE1_RESOLVENT_STEP["interior"], float(GARCIA_2021_CASE1_ROUNDING[:-1].max())
)
OUTER_SURFACE_TOLERANCE = tolerance_for(
    GARCIA_CASE1_RESOLVENT_STEP["surface"], float(GARCIA_2021_CASE1_ROUNDING[-1])
)


def solve_case1(n_r: int, n_mu: int):
    """Garcia 2021 Case 1 by the trajectory resolvent at ``(n_r, n_mu)``."""
    return solve_greens_function_sphere_mr_fixed_source(
        radii=GARCIA_2021_CASE1_RADII,
        sigma_t=np.array([[1.0], [0.5], [2.0]]),
        sigma_s=np.array([[[0.99]], [[0.3]], [[1.9]]]),
        external_source=np.array([[0.5], [1.0], [1.5]]),
        alpha=0.0,                 # vacuum BC
        n_r=n_r, n_mu=n_mu, n_traj_quad=64,
        max_iter=2000, tol=1e-7,
    )


def _phi_at(res, r: float) -> float:
    """The solved scalar flux at radius ``r``, from one cubic spline per region.

    At an interface radius (belonging to two regions) the two one-sided
    values are averaged. A radius outside [0, R] is refused.
    """
    outer = GARCIA_2021_CASE1_RADII
    if not 0.0 <= r <= outer[-1]:
        raise ValueError(f"r = {r} lies outside the sphere [0, {outer[-1]}]")
    pieces = _regionwise_cubic_spline(res.r_nodes, res.phi_g[0], res.region_at_node, len(outer))
    inner = np.concatenate([[0.0], outer[:-1]])
    values = [float(pieces[k](r)) for k, (a, b) in enumerate(zip(inner, outer)) if a <= r <= b]
    return float(np.mean(values))


@pytest.fixture(scope="module")
def garcia_2021_case1_solution():
    """Run the Variant α solver once for Garcia 2021 Case 1 and
    cache the result for all per-point assertions.
    """
    res = solve_case1(n_r=48, n_mu=24)
    return res


@pytest.mark.foundation
def test_garcia_case1_converged(garcia_2021_case1_solution):
    """Sanity: solver converges within iteration budget."""
    res = garcia_2021_case1_solution
    assert res.converged, (
        f"Garcia 2021 Case 1 did not converge in 2000 iter; "
        f"iter={res.iterations}"
    )


@pytest.mark.foundation
@pytest.mark.catches("ERR-090")
@pytest.mark.parametrize(
    "r, phi_garcia",
    [
        (r, p) for r, p in zip(GARCIA_2021_CASE1_R, GARCIA_2021_CASE1_PHI)
    ],
    ids=[f"r={r}" for r in GARCIA_2021_CASE1_R],
)
def test_garcia_case1_phi_matches_at_point(
    garcia_2021_case1_solution, r, phi_garcia,
):
    r"""Per-point flux cross-check against Garcia 2021 Table 5.

    The trajectory resolvent's flux (factor-of-2 convention applied) must
    match Garcia's converged ppP_N value within :data:`INTERIOR_TOLERANCE`
    inside the sphere and :data:`OUTER_SURFACE_TOLERANCE` at r = R; the
    derivation of both is in the module docstring.
    """
    res = garcia_2021_case1_solution
    phi_variant_alpha = _phi_at(res, r)
    phi_expected = GARCIA_TO_VARIANT_ALPHA_FACTOR * phi_garcia
    rel_err = abs(phi_variant_alpha - phi_expected) / phi_expected
    tol = OUTER_SURFACE_TOLERANCE if abs(r - GARCIA_2021_CASE1_RADII[-1]) < 1e-9 else INTERIOR_TOLERANCE
    print(f"r={r}: rel_err={rel_err:.2e} (tolerance {tol:.0e})")
    assert rel_err < tol, (
        f"Garcia Case 1 r={r}: Variant α = {phi_variant_alpha:.4f}, "
        f"expected = {phi_expected:.4f} (Garcia × 0.5 = "
        f"{phi_garcia} × 0.5), rel_err = {rel_err:.2e} > {tol:.0e}"
    )


@pytest.mark.foundation
def test_garcia_case1_qualitative_features(garcia_2021_case1_solution):
    r"""Qualitative shape properties of Garcia Case 1 solution.

    - Strong peaking in core region (most multiplication / source).
    - Through-region trends respect the sourcing pattern: φ decreases
      through core (high c, source 0.5), then increases-then-decreases
      through mid+outer regions (sources 1.0, 1.5).
    - φ vanishes at outer surface (vacuum BC).
    """
    res = garcia_2021_case1_solution
    phi = res.phi_g[0]

    # φ peaked near r=0.
    assert phi.argmax() < 5, (
        f"Garcia Case 1: φ should peak near r=0; got peak at "
        f"index {phi.argmax()}"
    )

    # φ at outer surface ≪ φ at centre (vacuum drains).
    centre_to_surface = phi[0] / phi[-1]
    assert centre_to_surface > 4.0, (
        f"Garcia Case 1: φ(centre)/φ(surface) = {centre_to_surface:.2f} "
        "should be > 4 (vacuum BC drains outer region)"
    )

    # φ positive everywhere.
    assert (phi > 0).all(), (
        f"Garcia Case 1: φ should be positive; min = {phi.min():.4e}"
    )
