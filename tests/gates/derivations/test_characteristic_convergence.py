"""The characteristic reference's self-convergence on its own resolution axes (the spec's D9).

P1 step (e1b) of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "Step (e): the audit
and four rulings"): the successors of the old trajectory-resolvent family's
convergence rows (its quadrature-floor, grazing-ray and off-diagonal-albedo rows,
19 CONV and 3 Q rows of the spec's roster), which step (e2) deletes with the
family. Verification spec ``scratch/characteristic_architecture/p1_verification_spec.md``,
row D9; the row-by-row map is ``scratch/characteristic_architecture/p1_step_e/ta_e1b/README.md``.

**The claim.** On each of the reference's resolution axes, varied alone from the
door's working point (degree 3, 2 grading layers of ratio 0.4, 8 points per
line-measure piece, 12 points per arc-length piece), the steps of k between
successive rungs contract strictly until they reach the solve's rounding, and
the working point's error, estimated per axis from its step up and the step below
(:func:`~tests.gates.derivations._ladder_rules.geometric_error`) and summed over
the axes, is below the old family's measured floor at the same kind of fixture
(the spec's acceptance: the old off-diagonal hollow sphere 3e-4 to 4e-4). A ladder
estimates an error and never certifies one (the user's step-5 ruling of
2026-10-03); the row's teeth are the contraction, which a defect that converges
to the wrong limit cannot break and a defect that breaks the rule's order does.

**The axes** (``TransportResolution(line_points, points, inner_points)``):
``degree`` (the panels' polynomial degree, with the projection points following
it), ``layers`` (the geometric grading toward every wall and interface), ``line``
(the points per piece of the line measure: the impact parameter, the polar angle,
the slab's cosine, which is where grazing lines live) and ``along`` (the points per
arc-length piece along each line and its inner rule).

**The fixtures**: two-group heterogeneous bodies PU2 | ABS | PU2 (fissile, a
downscattering absorber between) at intermediate albedos, the regime the old
off-diagonal rows posed: the solid sphere (0, 0.5, 1.5, 2.0) under 0.6, the
hollow sphere (0.4, 0.5, 1.5, 2.0) under (0.4, 0.7), the slab (0, 0.4, 1.5, 2.3)
under (0.4, 0.7); the cylinder (0, 1) under 0.6 and the annulus (0.4, 1.4) under
(0.4, 0.7), one region of PU2 each, on the angular axes at degree 2 (slow; the
three-region cylinder and annulus cost more than 10 min a solve, and their degree
axis is not ladder-gated here: a stated gap). ``[M]`` 2026-10-10 (``probes/probe_d9.log``): every axis contracts
on the sphere, hollow sphere and slab; the largest per-axis step at the working
point is 1.7e-6 (the layers, hollow sphere).
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.continuous.characteristic import CharacteristicDerivation, Resolution, TransportResolution
from orpheus.geometry import StructuredGeometry
from orpheus.geometry.boundary import AlbedoBoundary, SpecularReturn, VacuumInflow
from orpheus.numerics.observable import Eigenvalue
from orpheus.numerics.question import Eigen
from orpheus.numerics.traced_memo import bypass
from orpheus.specification.specification import GeometrySpecification
from tests.gates.derivations._ladder_rules import geometric_error
from tests.gates.derivations.test_characteristic_system import _ABS, _PU2

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_SYSTEM = "tests/gates/derivations/test_characteristic_system.py::"
_ALBEDOS = "tests/gates/derivations/test_characteristic_albedos.py::"
_D5 = _SYSTEM + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio"
_D8 = _ALBEDOS + "test_k_increases_strictly_with_each_walls_albedo_from_the_vacuum_body_toward_k_inf"

_K = CellCoefficient.every(Channel.FISSION_EMISSION)
_MR3, _MR3H, _SLB3 = (0.0, 0.5, 1.5, 2.0), (0.4, 0.5, 1.5, 2.0), (0.0, 0.4, 1.5, 2.3)


def _law(albedo: float):
    return VacuumInflow() if albedo == 0.0 else AlbedoBoundary(albedo, SpecularReturn(axis="x"))


_BODIES = {
    "sphere": lambda: StructuredGeometry.sphere(_MR3, (0, 1, 0), outer=_law(0.6)),
    "hollow": lambda: StructuredGeometry.sphere(_MR3H, (0, 1, 0), inner=_law(0.4), outer=_law(0.7)),
    "slab": lambda: StructuredGeometry.slab(_SLB3, (0, 1, 0), left=_law(0.4), right=_law(0.7)),
    # the cylinder and the annulus are one region (PU2): ``[M]`` 2026-10-10, a single three-region two-group annulus
    # solve at degree 2 ran past 10 min, a one-region one 10 s (``probes/probe_d8b.log``)
    "cylinder": lambda: StructuredGeometry.cylinder((0.0, 1.0), (0,), outer=_law(0.6)),
    "annulus": lambda: StructuredGeometry.cylinder((0.4, 1.4), (0,), inner=_law(0.4), outer=_law(0.7)),
}


def _resolution(degree: int = 3, layers: int = 2, line: int = 8, along: int = 12) -> Resolution:
    return Resolution(degree, layers, 0.4, TransportResolution(line, along, along), 2 * (degree + 1))


#: axis -> (the rungs, the working point's index among them). The working point is the door's default.
_AXES = {
    "degree": ([dict(degree=p) for p in (2, 3, 4, 5, 6)], 1),
    "layers": ([dict(layers=n) for n in (0, 1, 2, 3, 4)], 2),
    "line": ([dict(line=n) for n in (4, 6, 8, 12, 16)], 2),
    "along": ([dict(along=n) for n in (6, 8, 12, 16, 24)], 2),
}
#: The annulus's ladders at degree 2 (its rung-3 solve costs about 190 s): the angular axes only.
_ANNULUS_AXES = {
    "line": ([dict(degree=2, line=n) for n in (4, 6, 8, 12)], 2),
    "along": ([dict(degree=2, along=n) for n in (6, 8, 12, 16)], 2),
}
#: The cylinder's ladders here, as the annulus's: the angular axes at degree 2 on one region. ``[M]`` 2026-10-10: on a
#: three-region two-group cylinder the degree axis to 4 and the angular axis at degree 3 each ran past 30 min. The
#: degree axis on a multi-region cylinder is its own row,
#: ``test_k_contracts_in_the_panel_degree_on_a_two_region_cylinder`` (D5's flat closed bodies cannot see the degree
#: axis, and ``_characteristic_ladders.ABA_CYLINDER_K`` is a table whose only live gate was the corroboration file).
_CYLINDER_AXES = _ANNULUS_AXES
#: A step below this, relative, is the solve's rounding (the dense pencil's eigenvalue on a few hundred unknowns):
#: past it the ladder has converged and the steps need not contract. ``[M]``: the converged steps sit at 2e-16 to
#: 1.4e-14 on every axis of the three cheap bodies.
_ROUNDING = 1e-12
#: The old family's measured floor at an off-diagonal albedo (``test_off_diagonal_intermediate_alpha_hollow_sphere``:
#: 3e-4 to 4e-4), the spec's D9 acceptance: the estimated error at the working point is at most this.
_ACCEPTANCE = 3e-4


#: The angular axis (``line``) converges exponentially: each step above the rounding is at most this fraction of
#: the step before it. ``[M]`` 2026-10-10 (``probes/probe_d9.log``): 1.7e-3 (sphere), 1.4e-3 (hollow sphere),
#: 2.8e-4 (slab); with the tangency substitution dropped (the battery's ``tangency-substitution-dropped``) the
#: ratios are 0.29, 0.57, 0.26 (``probes/probe_d9_arm.py``): algebraic. The other axes contract more slowly
#: (degree 0.08 to 0.31, layers 0.01 to 0.07) and carry the contraction leg only.
_ANGULAR_RATE = 1e-2


@lru_cache(maxsize=None)
def _k(body: str, resolution: Resolution) -> float:
    spec = GeometrySpecification(Materials({0: _PU2, 1: _ABS}), _BODIES[body](), Eigen(_K))
    with bypass():
        return float(CharacteristicDerivation(spec, resolution).evaluate(Eigenvalue()).value)


def _contraction(body: str, rungs: list[dict]) -> tuple[list[float], np.ndarray]:
    ks = [_k(body, _resolution(**r)) for r in rungs]
    return ks, np.abs(np.diff(ks)) / abs(ks[-1])


def _assert_contracts(steps: np.ndarray, what: str) -> None:
    """Each step above the rounding is strictly below the one before it; once a step reaches the rounding, every
    later one stays there."""
    for i in range(1, len(steps)):
        if steps[i - 1] < _ROUNDING:
            assert steps[i] < _ROUNDING, (what, steps)
        else:
            assert steps[i] < steps[i - 1], (what, steps)


def _working_point_error(steps: np.ndarray, working: int) -> float:
    """The working point's error on one axis: its step up over one minus the ratio to the step below
    (``geometric_error``); the rounding where the step up is already there."""
    up = steps[working]
    if up < _ROUNDING:
        return _ROUNDING
    return geometric_error(up, steps[working - 1])


_CASES = [
    pytest.param("sphere", _AXES, id="sphere"),
    pytest.param("hollow", _AXES, id="hollow-sphere"),
    pytest.param("slab", _AXES, id="slab", marks=pytest.mark.slow),
    pytest.param("cylinder", _CYLINDER_AXES, id="cylinder", marks=pytest.mark.slow),
    pytest.param("annulus", _ANNULUS_AXES, id="annulus", marks=pytest.mark.slow),
]


@pytest.mark.l2
@pytest.mark.verifies("peierls-greens-cylinder-mr-quadrature-convergence")
@pytest.mark.parametrize(("body", "axes"), _CASES)
@pytest.mark.rests_on(_D5, _D8)
def test_k_contracts_on_every_resolution_axis_and_the_working_point_is_within_the_old_floor(body: str, axes) -> None:
    """[D9] On each axis, varied alone from the working point, the steps of k contract strictly down to the
    solve's rounding, the angular axis by at least 100x a step, and the working point's error summed over the
    axes is below the old family's floor.

    It succeeds the old family's convergence rows (the quadrature floors,
    the grazing-ray stability rows, the quadrature-order ladders and the
    off-diagonal intermediate-albedo rows of the sphere, hollow sphere, slab,
    asymmetric slab, cylinder and annulus), each of which asserted a floor or
    a finite k on its own ladder. A grazing-ray instability shows as a
    non-contracting ``line`` axis. First reds (``[M]`` 2026-10-10, battery
    ``scratch/characteristic_architecture/p1_step_e/ta_e1b/battery``): the
    impact rule's substitution y = sqrt(b^2 - r_k^2) dropped, so the angular
    axis converges algebraically (arm ``tangency-substitution-dropped``, the
    rate leg); the panels' grading ignored, so the degree axis stops
    contracting (arm ``layers-ignored``). Declared blind, measured green: the
    grading toward the tangencies dropped while the substitution is kept
    (arm ``tangency-ungraded``): the substitution alone carries the rate.
    """
    estimate = 0.0
    for axis, (rungs, working) in axes.items():
        ks, steps = _contraction(body, rungs)
        assert np.all(np.isfinite(ks)), (axis, ks)
        _assert_contracts(steps, f"{body}, {axis}")
        if axis == "line":
            for before, after in zip(steps, steps[1:], strict=False):
                if after >= _ROUNDING:
                    assert after < _ANGULAR_RATE * before, (body, "the angular axis is not exponential", steps)
        estimate += _working_point_error(steps, working)
    assert estimate < _ACCEPTANCE, estimate



def _two_region_cylinder_k(degree: int) -> float:
    """k of the vacuum cylinder of fuel A | moderator B at one group, (0, 0.6, 1.0) cm, at panel degree ``degree``."""
    from tests.gates.sn.verification.analytical._aba_reference import isotropic_mixture

    spec = GeometrySpecification(Materials({0: isotropic_mixture("A", "1g"), 1: isotropic_mixture("B", "1g")}),
                                 StructuredGeometry.cylinder((0.0, 0.6, 1.0), (0, 1), outer=_law(0.0)), Eigen(_K))
    resolution = Resolution(degree, 1, 0.4, TransportResolution(8, 8, 8), 2 * (degree + 1))
    with bypass():
        return float(CharacteristicDerivation(spec, resolution).evaluate(Eigenvalue()).value)


#: ``[M]`` 2026-10-10 (``probes/probe_cyl_cost.py``): degree 2 and 3 at 41 s and 48 s, k 0.58054060, 0.58054171.
_CYLINDER_DEGREES = (2, 3, 4, 5)


@pytest.mark.l2
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-cylinder-mr-quadrature-convergence")
@pytest.mark.rests_on(_D5)
def test_k_contracts_in_the_panel_degree_on_a_two_region_cylinder() -> None:
    """[D9, the degree axis on a multi-region cylinder] Fuel A | moderator B (0, 0.6, 1.0), one group, vacuum: k at
    degrees 2 to 5 (one grading layer, 8 points per piece) has strictly contracting steps down to the rounding, and
    the degree-3 error estimated from them (``geometric_error``) is below the old floor 3e-4.

    qa's ruling of 2026-10-10: a live gate on the multi-region cylinder's
    degree axis, which D5's flat closed bodies cannot see (a flat emission is
    in every degree's span) and whose step-(d) table ``ABA_CYLINDER_K`` loses
    its only live gate in (e2). The fixture is one group and two regions so
    each solve takes about a minute; the interface makes the emission's
    derivative jump, which the per-region panels must resolve. ``[M]``
    2026-10-10: 206 s. First red by value: the Gram matrix scaled by
    1 + 1e-6 (p + 1), an error growing with the degree (arm
    ``mass-degree-drift``: the steps stop contracting). Measured green: the
    grading ignored (``layers-ignored``; one layer at ratio 0.4 or none both
    converge in degree here); an under-integrated Gram (``mass-two-points``)
    reds structurally (a singular mass, ``NoFundamentalMode``).
    """
    ks = [_two_region_cylinder_k(p) for p in _CYLINDER_DEGREES]
    steps = np.abs(np.diff(ks)) / abs(ks[-1])
    _assert_contracts(steps, "two-region cylinder, degree")
    assert _working_point_error(steps, 1) < _ACCEPTANCE, (ks, steps)
