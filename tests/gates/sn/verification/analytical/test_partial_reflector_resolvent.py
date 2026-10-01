r"""ERR-094's L1 rung: the SN eigenvalue of a body with a partially reflecting
face, against the trajectory resolvent.

A face returning the fraction :math:`0 < \alpha < 1` of its outflow
specularly leaks :math:`(1 - \alpha) J^+`. Until 2026-09-30 the SN eigenvalue
counted the leakage only of faces returning nothing, so every such body
reported a k as if the face were a mirror (ERR-094). The estimator's identity
with the posed problem is gated in
``tests/gates/sn/eigenvalue/test_keff_estimator_gate.py::TestPartialReturnLeakage``;
those rows share the realized boundary operator with the solver. The rows here
compare the eigenvalue itself with a reference that shares nothing with SN
above the trusted-library line.

The reference is the trajectory resolvent, the semi-analytical pillar
(``vv-principles``, the three pillars): the integral transport equation along
characteristics, with the specular reflectivity entering a boundary-to-boundary
resolvent (rank 2 on the slab, one surface on the sphere). Its corners
(:math:`\alpha \in \{0, 1\}`) and the method of images are verified in
``tests/gates/derivations/test_peierls_greens_function_slab_asymmetric_solver.py``;
at intermediate :math:`\alpha` the evidence is the convergence of the two
methods to one value, measured below.

Tolerance. Each row compares an SN fixture with a resolvent fixture, and each
carries a discretization error. The band is the sum of the two errors, each
measured against the joint limit of both ladders, rounded up; it is not the
observed agreement. ``[M]`` 2026-09-30, two groups (mixture A,
:math:`k_\infty = 1.875`), SN with ``keff_tol=1e-10``:

=========  ==================================  ==================================
row        SN ladder (cells / Gauss-Legendre)  resolvent ladder
=========  ==================================  ==================================
slab       20/8 0.8303641, 40/16 0.8336481,    (16, 24, 32) 0.8343101,
L = 4 cm   80/16 0.8337041, 80/32 0.8344518,   x2 0.8346122, x3 0.8346710,
0.3 | 0.7  160/32 0.8344658, 320/64 0.8346559  x4 0.8346920 (n_x, n_mu, n_traj)
sphere     20/8 0.8831685, 40/16 0.8825984,    (12, 12, 32) 0.8829013,
R = 4 cm   80/32 0.8822111, 160/64 0.8820750   (24, 24, 64) 0.8821487,
0.7                                            (36, 36, 96) 0.8820581 (n_r, n_mu, n_traj)
=========  ==================================  ==================================

The finest members agree to 4.3e-5 (slab) and 1.9e-5 (sphere), relative. At
the fixtures used (slab: SN 80/32, resolvent x1; sphere: SN 40/16, resolvent
(24, 24, 64)) the errors against the joint limit are 3.0e-4 + 4.6e-4 (slab)
and 6.0e-4 + 0.9e-4 (sphere): the band is 1e-3 for both. The defect is three
orders above it: with the leakage predicate reverted, SN reads 1.82 on both
bodies (+118 % and +107 %).

What the band cannot see: the curvilinear corner of ERR-094 (the
off-quadrature :math:`\mu = \pm 1` ray re-emitting its full outflow) moves
the sphere's k by an amount of the order of the angular discretization
error, which vanishes as the quadrature is refined (``[M]`` 2026-09-30, the
archivist's ladder on a one-group sphere of radius 2 cm with
:math:`\alpha = 0.7`: scaled and unscaled corners 2.3e-4 apart, relative, at
GL16 and 80 cells, 1.2e-6 at GL64 and 320 cells). Measured on this module
(``scratch/boundary_ontology/battery_err094.md``, 2026-09-30): the sphere row
stays green under the unscaled corner, and reds only when the corner is
silenced altogether (a control arm). The corner's catchers are the
operator-level gate
``tests/gates/sn/operators/test_psi_half_coupling.py::TestB_b_RayBoundary::test_partial_specular_corner_is_alpha_times_the_mirror``
and the :math:`\alpha = 0` edge of ``TestPartialReturnLeakage``.

What these rows see that the estimator rows cannot: an operator realizing the
wrong amplitude. With the SN realizer returning :math:`\alpha^2` instead of
:math:`\alpha` both rows red while every map-ratio row stays green (the same
battery).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.continuous.trajectory_resolvent.greens_function import (
    solve_greens_function_sphere_mg,
)
from orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab_asymmetric import (
    solve_greens_function_slab_asymmetric_mg,
)
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.boundary import AlbedoBoundary, ReflectiveBoundary, SpecularReturn
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn

pytestmark = pytest.mark.l1

_ESTIMATOR = (
    "tests/gates/sn/eigenvalue/test_keff_estimator_gate.py::TestPartialReturnLeakage::"
    "test_reported_k_is_the_posed_eigenvalue"
)
_RESOLVENT_IMAGES = (
    "tests/gates/derivations/test_peierls_greens_function_slab_asymmetric_solver.py::"
    "test_method_of_images_reflective_vacuum_equals_double_vacuum"
)
_MIRROR_ACTION = (
    "tests/gates/geometry/test_reemission_closure.py::"
    "TestSpecularAgainstAnIndependentExpression::"
    "test_matches_the_hand_written_mirror_gather"
)

#: The band: the two fixtures' discretization errors against the joint
#: limit of both ladders, summed and rounded up (module docstring).
_RTOL = 1e-3

_FUEL = get_mixture("A", "2g")


def _cross_sections():
    """Mixture A as the resolvent reads it: the P0 transfer ``[g_from, g_to]``
    (SN is solved at scattering order 0, so both read the same kernel)."""
    return (
        np.asarray(_FUEL.SigT, dtype=float), _FUEL.SigS[0].toarray(),
        np.asarray(_FUEL.SigP, dtype=float), np.asarray(_FUEL.chi, dtype=float),
    )


def _sn_k(mesh, n_ordinates: int) -> float:
    solution = solve_sn(
        {0: _FUEL}, mesh, Quadrature.gauss_legendre(n_ordinates),
        keff_tol=1e-10, flux_tol=1e-9, inner_tol=1e-12, max_inner=4000,
    )
    return float(solution.outcome.keff)


@pytest.mark.rests_on(_ESTIMATOR, _RESOLVENT_IMAGES, _MIRROR_ACTION)
@pytest.mark.verifies("sn-keff-update")
@pytest.mark.verifies("sn-leakage-functional")
@pytest.mark.catches("ERR-094")
def test_slab_with_two_partial_reflectors():
    r"""A homogeneous slab of 4 cm with ``ReflectiveBoundary("x", 0.3)`` on
    the left and ``AlbedoBoundary(0.7, SpecularReturn("x"))`` on the right:
    two amplitudes and both spellings on one body. Claim kind: REFERENCE.

    First red, ``[M]`` 2026-09-30 with the leakage predicate reverted to
    ``response_kernel.is_zero``: SN 1.82163 against the resolvent 0.83431
    (+118 %).
    """
    sigma_t, sigma_s, nu_sigma_f, chi = _cross_sections()
    reference = solve_greens_function_slab_asymmetric_mg(
        4.0, sigma_t, sigma_s, nu_sigma_f, chi,
        alpha_left=0.3, alpha_right=0.7,
        n_x=16, n_mu=24, n_traj_quad=32, max_iter=500, tol=1e-12,
    )
    if not reference.converged:
        pytest.fail("the resolvent did not converge: no reference to compare")
    mesh = Mesher(StructuredGeometry.slab(
        (0.0, 4.0), (0,),
        left=ReflectiveBoundary("x", 0.3),
        right=AlbedoBoundary(0.7, SpecularReturn("x")),
    )).partition(CellsByCount.uniform_width(80)).mesh
    k = _sn_k(mesh, 32)
    np.testing.assert_allclose(
        k, reference.k_eff, rtol=_RTOL,
        err_msg="SN's k of a slab with partially reflecting faces is not the "
                "resolvent's within the two methods' discretization bands "
                "(ERR-094: the (1 - alpha) J+ of each face is leakage).",
    )


@pytest.mark.rests_on(_ESTIMATOR, _MIRROR_ACTION)
@pytest.mark.verifies("sn-keff-update")
@pytest.mark.verifies("sn-leakage-functional")
@pytest.mark.catches("ERR-094")
def test_sphere_with_a_partial_reflector():
    r"""A homogeneous sphere of radius 4 cm with ``ReflectiveBoundary("x",
    0.7)`` on its surface. Claim kind: REFERENCE.

    First red, ``[M]`` 2026-09-30 with the leakage predicate reverted to
    ``response_kernel.is_zero``: SN about 1.82 against the resolvent 0.88215
    (+107 %).
    """
    sigma_t, sigma_s, nu_sigma_f, chi = _cross_sections()
    reference = solve_greens_function_sphere_mg(
        4.0, sigma_t, sigma_s, nu_sigma_f, chi,
        alpha=0.7, n_r=24, n_mu=24, n_traj_quad=64, max_iter=500, tol=1e-11,
    )
    if not reference.converged:
        pytest.fail("the resolvent did not converge: no reference to compare")
    mesh = Mesher(StructuredGeometry.sphere(
        (0.0, 4.0), (0,), outer=ReflectiveBoundary("x", 0.7),
    )).partition(CellsByCount.uniform_width(40)).mesh
    k = _sn_k(mesh, 16)
    np.testing.assert_allclose(
        k, reference.k_eff, rtol=_RTOL,
        err_msg="SN's k of a sphere with a partially reflecting surface is not "
                "the resolvent's within the two methods' discretization bands "
                "(ERR-094).",
    )
