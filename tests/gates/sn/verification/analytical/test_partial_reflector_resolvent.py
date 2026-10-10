r"""ERR-094's L1 rung: the SN eigenvalue of a body with a partially reflecting
face, against the characteristic reference.

A face returning the fraction :math:`0 < \alpha < 1` of its outflow
specularly leaks :math:`(1 - \alpha) J^+`. Until 2026-09-30 the SN eigenvalue
counted the leakage only of faces returning nothing, so every such body
reported a k as if the face were a mirror (ERR-094). The estimator's identity
with the posed problem is gated in
``tests/gates/sn/eigenvalue/test_keff_estimator_gate.py::TestPartialReturnLeakage``;
those rows share the realized boundary operator with the solver. The rows here
compare the eigenvalue itself with a reference whose solve shares nothing
with SN above the trusted-library line. One datum is shared: both sides read
the face's response factor from the same object, the law's
``SpecularReturn.kernel`` (the reference through ``law.response_kernel`` in
``characteristic/walls.py``, SN through its leakage predicate and its
corner). It is part of the posed problem, not an independent check: ``[M]``
2026-10-10 (qa, ``p1_step_d/qa_d/half_alpha_plugin.py``), alpha := alpha / 2
there moves both sides' k together and the slab row stays green. Its
external pin is ``tests/gates/derivations/test_characteristic_walls.py``,
red 3 times under that mutation.

The reference is the characteristic reference (the semi-analytical pillar,
``vv-principles``, the three pillars): the transport equation integrated along
the lines of the body, Galerkin over lines, the partial mirror entering the
boundary resolvent of each line (``test_characteristic_closure.py``, the
unfolded wall-by-wall sum, is the method of images it rests on). It is posed
from the same specification SN solves, mixture A at P0
(``_aba_reference.isotropic_mixture``): SN solves at scattering order 0, so the
library mixture's P1 moment, which the reference refuses, is read by neither
side (asserted in each row). Until P1 step (d) of
``.claude/plans/characteristic_reference_architecture.md`` (2026-10-09) it was
the trajectory resolvent (rank 2 on the slab, one surface on the sphere).

Tolerance. Each row compares an SN fixture with a reference fixture, and each
carries a discretization error. The band is the sum of the two errors, each
measured against a limit, rounded up (``_ladder_rules.summed_tolerance``,
the user's ruling 2 of 2026-10-09: these rows keep their sum rule, which has
no factor 2): it is not the observed agreement. ``[M]`` 2026-09-30, two groups
(mixture A, :math:`k_\infty = 1.875`), SN with ``keff_tol=1e-10``, and the
retired resolvent's ladders:

=========  ==================================  ==================================
row        SN ladder (cells / Gauss-Legendre)  resolvent ladder (retired)
=========  ==================================  ==================================
slab       20/8 0.8303641, 40/16 0.8336481,    (16, 24, 32) 0.8343101,
L = 4 cm   80/16 0.8337041, 80/32 0.8344518,   x2 0.8346122, x3 0.8346710,
0.3 | 0.7  160/32 0.8344658, 320/64 0.8346559  x4 0.8346920 (n_x, n_mu, n_traj)
sphere     20/8 0.8831685, 40/16 0.8825984,    (12, 12, 32) 0.8829013,
R = 4 cm   80/32 0.8822111, 160/64 0.8820750   (24, 24, 64) 0.8821487,
0.7                                            (36, 36, 96) 0.8820581 (n_r, n_mu, n_traj)
=========  ==================================  ==================================

At the SN fixtures used (slab 80/32, sphere 40/16) SN's error against the
joint limit of both ladders is 3.0e-4 (slab) and 6.0e-4 (sphere): the
:data:`_SN_ERROR` terms. The reference's term is its own ladder's estimate at
p = 5 (``_characteristic_ladders.reference_error``): about 1e-9 (slab) and
6e-8 (sphere), three orders below the retired resolvent's 4.6e-4 and 0.9e-4.
The bands are therefore 4e-4 (slab) and 7e-4 (sphere), tightened from 1e-3.
``[M]`` 2026-10-09 the readings are 3.2e-4 (slab) and 6.5e-4 (sphere): SN's
error against the characteristic reference, 7 % above the 6.0e-4 measured
against the old joint limit on the sphere, so the sphere row holds a margin of
1.07 (the ruling records it). The defect is three orders above either band:
with the leakage predicate reverted, SN reads 1.82 on both bodies (+118 % and
+107 %).

What the band cannot see: the curvilinear corner of ERR-094 (the
off-quadrature :math:`\mu = \pm 1` ray re-emitting its full outflow) moves
the sphere's k by an amount of the order of the angular discretization
error, which vanishes as the quadrature is refined (``[M]`` 2026-09-30, the
archivist's ladder on a one-group sphere of radius 2 cm with
:math:`\alpha = 0.7`: scaled and unscaled corners 2.3e-4 apart, relative, at
GL16 and 80 cells, 1.2e-6 at GL64 and 320 cells). Measured on this module
(``scratch/boundary_ontology/battery_err094.md``, 2026-09-30, against the
retired resolvent at the band 1e-3): the sphere row stays green under the
unscaled corner, and reds only when the corner is silenced altogether (a
control arm). At step (d)'s 7e-4 band, ``[M]`` 2026-10-10 (qa,
``p1_step_d/qa_d/corner_*.log``): re-dropping the unscaled corner (alpha := 1
in ``_reflect_corner``, 1465 activations) moves SN's k by 9.4e-4, but TOWARD
the reference: the reading goes from 6.54e-4 to 2.86e-4. The band cannot see
the corner at this fixture. The corner's catchers are the
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

import inspect

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.mesh import CellsByCount, Mesher
from orpheus.numerics.observable import Eigenvalue
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import solve_sn
from tests.gates.derivations._characteristic_ladders import (
    WORKING_DEGREE,
    partial_reflector_specification,
    reference_error,
    rung,
)
from tests.gates.derivations._ladder_rules import summed_tolerance
from tests.gates.sn.verification.analytical._aba_reference import isotropic_mixture

pytestmark = pytest.mark.l1

_ESTIMATOR = (
    "tests/gates/sn/eigenvalue/test_keff_estimator_gate.py::TestPartialReturnLeakage::"
    "test_reported_k_is_the_posed_eigenvalue"
)
#: The reference's method of images: a partial mirror's cycle is its unfolded wall-by-wall sum, at two distinct
#: amplitudes (the slab's pair) on a rank-2 line. Until step (d) the retired resolvent's
#: ``test_peierls_greens_function_slab_asymmetric_solver.py::test_method_of_images_reflective_vacuum_equals_double_vacuum``.
_REFERENCE_IMAGES = (
    "tests/gates/derivations/test_characteristic_closure.py::"
    "test_the_inflow_is_the_unfolded_wall_by_wall_sum[a0.3_0.6-2]"
)
_MIRROR_ACTION = (
    "tests/gates/geometry/test_reemission_closure.py::"
    "TestSpecularAgainstAnIndependentExpression::"
    "test_matches_the_hand_written_mirror_gather"
)

#: SN's relative error in k at each row's fixture (slab 80 cells / GL32, sphere 40 / GL16), against the joint
#: limit of the SN and the retired resolvent ladders: ``[M]`` 2026-09-30, the module docstring's table.
_SN_ERROR = {"slab": 3.0e-4, "sphere": 6.0e-4}


def _band(body: str) -> float:
    """The row's relative band: SN's error plus the reference's at its working point, rounded up (module docstring)."""
    return summed_tolerance(_SN_ERROR[body], reference_error(f"partial_reflector_{body}"))


_FUEL = get_mixture("A", "2g")


def _assert_one_problem() -> None:
    """SN and the reference pose one problem: the reference's mixture is SN's at P0, bit for bit, and SN solves at
    scattering order 0, its default (``_sn_k`` passes none), so the P1 moment SN is handed is read by neither side."""
    posed = isotropic_mixture("A")
    np.testing.assert_array_equal(np.asarray(posed.SigT), np.asarray(_FUEL.SigT))
    np.testing.assert_array_equal(posed.SigS[0].toarray(), _FUEL.SigS[0].toarray())
    np.testing.assert_array_equal(np.asarray(posed.SigP), np.asarray(_FUEL.SigP))
    np.testing.assert_array_equal(np.asarray(posed.chi), np.asarray(_FUEL.chi))
    assert len(posed.SigS) == 1 and len(_FUEL.SigS) == 2, "the premise: the library mixture carries a P1 moment"
    assert inspect.signature(solve_sn).parameters["scattering_order"].default == 0


def _reference_k(body: str) -> float:
    """The characteristic reference's k on the row's body at its working point."""
    from orpheus.derivations.continuous.characteristic import characteristic_reference

    resolution = rung(WORKING_DEGREE[f"partial_reflector_{body}"])
    return float(characteristic_reference(partial_reflector_specification(body), resolution).read(Eigenvalue()).value)


def _sn_k(mesh, n_ordinates: int) -> float:
    solution = solve_sn(
        {0: _FUEL}, mesh, Quadrature.gauss_legendre(n_ordinates),
        keff_tol=1e-10, flux_tol=1e-9, inner_tol=1e-12, max_inner=4000,
    )
    return float(solution.outcome.keff)


@pytest.mark.rests_on(_ESTIMATOR, _REFERENCE_IMAGES, _MIRROR_ACTION)
@pytest.mark.verifies("sn-keff-update")
@pytest.mark.verifies("sn-leakage-functional")
@pytest.mark.catches("ERR-094")
def test_slab_with_two_partial_reflectors():
    r"""A homogeneous slab of 4 cm with the specular wall
    ``AlbedoBoundary(0.3, SpecularReturn("x"))`` on the left and
    ``AlbedoBoundary(0.7, SpecularReturn("x"))`` on the right: two amplitudes
    on one body. Claim kind: REFERENCE. Until the reflective cleanup the left
    wall was spelled ``ReflectiveBoundary("x", 0.3)``, the same matrix
    (``[M]`` 2026-10-01, the carry fixture ``resolvent_slab``: SN k bitwise
    equal).

    Band 4e-4 relative (module docstring), reading 3.2e-4 ``[M]`` 2026-10-09.
    First red, ``[M]`` 2026-09-30 with the leakage predicate reverted to
    ``response_kernel.is_zero``: SN 1.82163 against the resolvent 0.83431
    (+118 %).
    """
    _assert_one_problem()
    specification = partial_reflector_specification("slab")
    mesh = Mesher(specification.geometry).partition(CellsByCount.uniform_width(80)).mesh
    k = _sn_k(mesh, 32)
    np.testing.assert_allclose(
        k, _reference_k("slab"), rtol=_band("slab"),
        err_msg="SN's k of a slab with partially reflecting faces is not the "
                "reference's within the two methods' discretization bands "
                "(ERR-094: the (1 - alpha) J+ of each face is leakage).",
    )


@pytest.mark.rests_on(_ESTIMATOR, _MIRROR_ACTION)
@pytest.mark.verifies("sn-keff-update")
@pytest.mark.verifies("sn-leakage-functional")
@pytest.mark.catches("ERR-094")
def test_sphere_with_a_partial_reflector():
    r"""A homogeneous sphere of radius 4 cm with the specular wall
    ``AlbedoBoundary(0.7, SpecularReturn("x"))`` on its surface (spelled
    ``ReflectiveBoundary("x", 0.7)`` until the reflective cleanup; ``[M]``
    2026-10-01, the carry fixture ``resolvent_sphere``: SN k bitwise equal).
    Claim kind: REFERENCE.

    Band 7e-4 relative (module docstring), reading 6.5e-4 ``[M]`` 2026-10-09: a margin of 1.07, because SN's own
    error at 40 cells and GL16 against this reference (6.5e-4) is 7 % above the 6.0e-4 the band's SN term carries,
    which was measured against the joint limit of the SN and retired-resolvent ladders (the user's ruling 2 of
    2026-10-09). First red, ``[M]`` 2026-09-30 with the leakage predicate reverted to
    ``response_kernel.is_zero``: SN about 1.82 against the resolvent 0.88215
    (+107 %).
    """
    _assert_one_problem()
    specification = partial_reflector_specification("sphere")
    mesh = Mesher(specification.geometry).partition(CellsByCount.uniform_width(40)).mesh
    k = _sn_k(mesh, 16)
    np.testing.assert_allclose(
        k, _reference_k("sphere"), rtol=_band("sphere"),
        err_msg="SN's k of a sphere with a partially reflecting surface is not "
                "the reference's within the two methods' discretization bands "
                "(ERR-094).",
    )
