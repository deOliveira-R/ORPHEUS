r"""The characteristic reference's error at each SN row's fixture, the SN side's residuals, and the command that re-measures each.

P1 step (d) of ``.claude/plans/characteristic_reference_architecture.md``
re-points the SN rows that compared against the trajectory resolvent onto the
characteristic reference (:func:`~orpheus.derivations.continuous.characteristic.characteristic_reference`).
Every tolerance of those rows is COMPUTED by the rules of
:mod:`tests.gates.derivations._ladder_rules` from the tables here, except one
term: the ERR-094 partial-reflector rows' SN error (``_SN_ERROR`` in
``test_partial_reflector_resolvent.py``, measured with that file's SN ladder);
their reference term is here. No tolerance is typed.

"The door default" is the door gates' resolution, ``_RES`` of
``tests/gates/derivations/test_characteristic_reference.py``:
``Resolution(3, 2, 0.4, TransportResolution(8, 12, 12), 8)``, :func:`rung` (3)
here. The rows:

- ``tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py`` (the
  homogeneous edge rows, the A|B|A sphere's k and shape, the A|B|A cylinder's
  k and shape, the RECORD);
- ``tests/gates/sn/verification/analytical/test_l1_standoff_slab_cylinder.py``
  and ``tests/gates/sn/sweep/curvilinear/test_unified_matvec_cylinder.py``
  (the A|B|A cylinder at folded 4x8, and the RECORD);
- ``tests/gates/sn/verification/analytical/test_partial_reflector_resolvent.py``
  (ERR-094's slab and sphere).

**What an estimate is here.** A ladder ESTIMATES a reference's error; it can
refute a reference, never certify one (the user's step-5 ruling of
2026-10-03, #405 P2). The characteristic reference reads ``Uncertified`` (#566),
so no number here is a bound, and the rows comparing against it are
``compare_uncertified`` comparisons or strict xfails on the verbs' refusal.

**The ladder.** One joint ladder, every axis that enters an eigenvalue moving
together (panel degree, grading layers, the transport block's line and
traversal points, the projection's source points): :data:`RUNGS`, keyed by the
degree. ``source_points`` projects a symbolic source, detector or weight and
does not enter k. Each fixture is read at its WORKING POINT (the rung the SN
rows use) and at the rung ABOVE it, and the estimate is the working point's
error, from the step up:

- **geometric** (:func:`~tests.gates.derivations._ladder_rules.geometric_error`)
  where a third rung below the working point gives the contraction ratio of the
  steps: the hp ladder converges exponentially, so the steps shrink
  geometrically (``[M]`` ratios 0.06 to 0.30 in k), and the working point's
  error is its step up over one minus that ratio. Not a last step: the last
  step understates the error by that factor. On the p = 5 fixtures the
  ladder runs to p = 8 and each estimate covers the distance to p = 8
  (``test_crosscheck_harness.py::test_each_reference_estimate_covers_the_finest_measured_rung``).
  Two variants, each measured: the A|B|A sphere's SHAPE converges in pairs of
  degrees, so its steps are two rungs long (:data:`ABA_SPHERE_SHAPE_STEP`);
  the A|B|A cylinder climbs its PANEL ladder only, its transport axis added
  as its own measured step (:func:`cylinder_panel_rung`);
- **closed form** on the homogeneous edge bodies, whose k is :math:`k_\infty`
  of the medium: the distance to it is the error itself, and the estimate is
  the larger of that distance and the step up.

Re-measure with::

    python -O -m tests.gates.derivations._characteristic_ladders <ladder>

with ``<ladder>`` one of ``aba-sphere``, ``aba-cylinder [degree ...]``,
``partial-reflector``, ``edge``, ``sphere-sn``, ``cylinder-sn``, ``records`` (the old family's two ladders,
``old-sphere-reference`` and ``old-cylinder-reference``, were deleted with it in step (e2)). Run from the repository
root, serially (the host is shared). The reference rungs are read in process
(``traced_memo.bypass()``); the ``records`` command goes through the memo, as
the rows do, and leaves the entry they read.
"""

from __future__ import annotations

import sys
import time
from typing import TYPE_CHECKING

from tests.gates.derivations._ladder_rules import geometric_error

if TYPE_CHECKING:
    from orpheus.derivations.continuous.characteristic import Resolution
    from orpheus.specification.specification import GeometrySpecification

# ── the ladder ────────────────────────────────────────────────────────────


def rung(degree: int) -> Resolution:
    """The joint ladder's rung of panel degree ``degree``: (degree, layers, ratio, transport, source points)."""
    from orpheus.derivations.continuous.characteristic import Resolution, TransportResolution

    layers, transport, source_points = _RUNG_AXES[degree]
    return Resolution(degree, layers, 0.4, TransportResolution(*transport), source_points)


#: degree -> (layers, (line points, traversal points, inner points), source points). Degree 3 is the door's gates'
#: resolution (``test_characteristic_reference.py``, ``_RES``); degrees 3 to 5 are the rungs step (c) measured on
#: (``p1_step_c/explorer/g12_sphere.py``); 2 and 6 extend them by one step each way.
_RUNG_AXES = {
    2: (1, (4, 8, 8), 6),
    3: (2, (8, 12, 12), 8),
    4: (3, (12, 16, 16), 10),
    5: (4, (16, 20, 20), 12),
    6: (5, (20, 24, 24), 14),
    7: (6, (24, 28, 28), 16),
    8: (7, (28, 32, 32), 18),
}

#: The working point of each fixture, by degree: the sphere and slab bodies at p = 5 (seconds), the A|B|A cylinder
#: and the homogeneous edge bodies at the door default p = 3 (the cylinder's p = 3 solve is about 21 min).
WORKING_DEGREE = {
    "aba_sphere": 5,
    "aba_cylinder": 3,
    "partial_reflector_slab": 5,
    "partial_reflector_sphere": 5,
    "edge_sphere_2g": 3,
    "edge_cylinder_1g": 3,
}


# ── the A|B|A sphere: k, and the 80-ratio shape ──────────────────────────
#: ``aba-sphere``: degree -> k. ``[M]`` 2026-10-09 (``p1_step_d/ta_d/ladders_p5.log``), in process, 0.7 s at p = 3 to
#: 202 s at p = 8. The steps up shrink geometrically: 1.48e-8, 9.30e-11, 2.75e-11, 1.59e-13, 3.17e-13 (the last two at
#: the solve's rounding).
ABA_SPHERE_K: dict[int, float] = {
    3: 1.3810544625961212,
    4: 1.3810544420949256,
    5: 1.3810544419664763,
    6: 1.3810544419284978,
    7: 1.381054441928278,
    8: 1.3810544419278405,
}
#: ``aba-sphere``: the shape metric between two rungs, ``max |r_p - r_q| / M_q`` over the 80 fission-gauged cell
#: averages on the SN rows' 40-cell mesh (``_aba_reference.shape_observables``), q the finer rung. ``[M]`` 2026-10-09.
#: The shape converges in PAIRS of degrees: the steps up are 2.31e-4, 9.30e-6, 9.07e-6, 4.41e-7, 4.72e-7, so a step of
#: one rung does not contract from p = 4 to p = 5 (ratio 0.975) and its geometric estimate would read 3.6e-4; the
#: estimate takes steps of two rungs (3 -> 5 -> 7, ratio 0.041). The pairing is on the degree axis alone (``[M]``
#: ``p1_step_d/ta_d/probe_shape_axes.log``: at p = 5 raising the layers moves the shape 4.9e-11, the line and traversal
#: points below 2e-14, the degree 9.07e-6) and is not the reading's (source points 12 -> 48 move it 1.3e-15,
#: ``probe_shape_sp.log``).
ABA_SPHERE_SHAPE_STEP: dict[tuple[int, int], float] = {
    (3, 5): 2.2123e-4,
    (5, 7): 9.0042e-6,
    (5, 8): 9.0324e-6,
}

# ── the A|B|A cylinder: k and the 80-ratio shape, on its own panel ladder ──
#: The cylinder's working point is the door default (p = 3, about 21 min a solve), and the joint ladder's rung above
#: it costs hours: on the homogeneous edge cylinder the joint step p = 3 -> 4 multiplies the cost 7.6 times (29 s ->
#: 220 s, ``ladders_cheap.log``), so about 2.7 h here. The cylinder's ladder therefore moves the PANEL axes only
#: (degree and grading layers together, the transport block held at the door default's (8, 12, 12)), and the
#: transport axis enters the estimate through its own measured step: the per-axis probes one rung below the default
#: (``p1_step_d/explorer/cyl_rung.log``) moved k 3.4e-8 for the degree and 6.9e-11 for the transport block, and on
#: the A|B|A sphere at p = 5 the degree carries the whole step (``probe_shape_axes.log``).
def cylinder_panel_rung(degree: int) -> Resolution:
    """The A|B|A cylinder's panel ladder: degree p, p - 1 grading layers, the door default's transport block."""
    from orpheus.derivations.continuous.characteristic import Resolution, TransportResolution

    return Resolution(degree, degree - 1, 0.4, TransportResolution(*_RUNG_AXES[3][1]), 2 * (degree + 1))


#: ``aba-cylinder``: panel degree -> k. ``[M]`` 2026-10-09/10, one solve each, serially, in process
#: (``p1_step_d/ta_d/ladder_cylinder.log``): p = 2 in 375 s, p = 3 (the door default) in 1259 s, p = 4 in 3399 s.
#: The steps up are 4.00e-7 and 2.44e-8 (ratio 0.061). The p = 3 value agrees with step (c)'s
#: 1.231729452890254 (``p1_step_c/iterated.log``) to one ulp.
ABA_CYLINDER_K: dict[int, float] = {
    2: 1.2317299451673347,
    3: 1.2317294528902538,
    4: 1.231729422811464,
}
#: ``aba-cylinder``: the shape metric between two panel rungs, as :data:`ABA_SPHERE_SHAPE_STEP`'s, on the cylinder
#: rows' 40-cell mesh (``p1_step_d/ta_d/ladder_cylinder_steps.log``): ratio 0.21.
ABA_CYLINDER_SHAPE_STEP: dict[tuple[int, int], float] = {
    (2, 3): 8.5496e-4,
    (3, 4): 1.7907e-4,
}
#: The transport block's step one rung below the door default, (6, 9, 9) -> (8, 12, 12), relative in k: ``[M]``
#: 2026-10-09 (``p1_step_d/explorer/cyl_rung.log``, 543 s): the error estimate of the coarser rung, so an
#: over-estimate of the default's transport error, added to the panel estimate.
ABA_CYLINDER_TRANSPORT_STEP = 6.906e-11

# ── the ERR-094 partial reflectors: k ─────────────────────────────────────
#: ``partial-reflector``: body -> degree -> k (slab L = 4 cm, albedos 0.3 | 0.7; sphere R = 4 cm, albedo 0.7; mixture
#: A at P0, two groups). ``[M]`` 2026-10-09 (``p1_step_d/ta_d/probe_joint_high.log``, reproduced bit for bit by
#: ``ladders_p5.log``). Steps up, slab: 2.48e-7, 4.31e-9, 8.81e-10, 3.26e-12, 4.06e-12; sphere: 1.55e-5, 1.02e-6,
#: 5.72e-8, 2.61e-9, 1.28e-10.
PARTIAL_REFLECTOR_K: dict[str, dict[int, float]] = {
    "slab": {
        3: 0.8347191534884694,
        4: 0.8347193605152384,
        5: 0.8347193641158219,
        6: 0.8347193648512348,
        7: 0.8347193648539515,
        8: 0.8347193648573381,
    },
    "sphere": {
        3: 0.8820071473279709,
        4: 0.8820207901846953,
        5: 0.8820216933851701,
        6: 0.8820217437908983,
        7: 0.8820217460941235,
        8: 0.8820217462073154,
    },
}

# ── the homogeneous edge bodies: k against k_inf ──────────────────────────
#: ``edge``: body -> degree -> k (mixture A at P0, R = 2 cm, reflective): the sphere at two groups, the cylinder at one.
#: ``[M]`` 2026-10-09 (``p1_step_d/ta_d/ladders_cheap.log``): 0.1 s and 0.3 s (sphere), 29 s and 220 s (cylinder).
EDGE_K: dict[str, dict[int, float]] = {
    "sphere_2g": {3: 1.875000000000014, 4: 1.875000000000019},
    "cylinder_1g": {3: 1.499999999999997, 4: 1.5000000000000016},
}
#: ``edge``: body -> the medium's k_inf (``kinf_homogeneous``, the closed form; ``edge_kinf``).
EDGE_KINF: dict[str, float] = {"sphere_2g": 1.8750000000000004, "cylinder_1g": 1.5}
#: ``edge``: each edge snapshot's own relative distance from its medium's k_inf, frozen with the snapshot: the SN side's
#: residual on the edge rows. ``[M]`` 2026-10-09 (``p1_step_d/ta_d/edge_snapshots.log``).
EDGE_SN_RESIDUAL: dict[str, float] = {
    "sphere_2g_homogeneous_dd_n20": 4.4604e-12,
    "cyl_1g_homogeneous_folded_4x8_dd_n20": 1.4803e-16,
    "cyl_1g_homogeneous_folded_2x4_dd_n20": 2.9606e-16,
}
#: Each edge snapshot's body.
_EDGE_SNAPSHOT_BODY = {
    "sphere_2g_homogeneous_dd_n20": "sphere_2g",
    "cyl_1g_homogeneous_folded_4x8_dd_n20": "cylinder_1g",
    "cyl_1g_homogeneous_folded_2x4_dd_n20": "cylinder_1g",
}


def _relative_step(k: dict[int, float], low: int, high: int) -> float:
    return abs(k[high] - k[low]) / abs(k[high])


def reference_error(fixture: str) -> float:
    """The characteristic reference's relative error ESTIMATE in k at the fixture's working point (module docstring)."""
    p = WORKING_DEGREE[fixture]
    match fixture:
        case "aba_sphere":
            k = ABA_SPHERE_K
        case "aba_cylinder":
            return _relative_step(ABA_CYLINDER_K, p, p + 1) / (1.0 - largest_measured_ratio()) + ABA_CYLINDER_TRANSPORT_STEP
        case "partial_reflector_slab" | "partial_reflector_sphere":
            k = PARTIAL_REFLECTOR_K[fixture.removeprefix("partial_reflector_")]
        case "edge_sphere_2g" | "edge_cylinder_1g":
            body = fixture.removeprefix("edge_")
            k = EDGE_K[body]
            return max(abs(k[p] - EDGE_KINF[body]) / EDGE_KINF[body], _relative_step(k, p, p + 1))
        case _:
            raise KeyError(fixture)
    return geometric_error(_relative_step(k, p, p + 1), _relative_step(k, p - 1, p))


def largest_measured_ratio() -> float:
    """The largest contraction ratio of the steps measured at a working point on the p = 5 ladders (k of the three
    bodies, the A|B|A sphere's shape over pairs): 0.30, the A|B|A sphere's k (``[M]`` 2026-10-09).

    The A|B|A cylinder's estimate uses it in place of its own ratio (qa F3, 2026-10-10): the cylinder's step below
    its working point, 2 -> 3, moves the grading layers from 1 to 2, and that step is dominated by the grading, so its
    ratio (0.061 in k) says little about the steps above. The largest ratio any ladder showed is the cautious one.
    """
    ratios = [_relative_step(k, 5, 6) / _relative_step(k, 4, 5)
              for k in (ABA_SPHERE_K, *PARTIAL_REFLECTOR_K.values())]
    ratios.append(ABA_SPHERE_SHAPE_STEP[(5, 7)] / ABA_SPHERE_SHAPE_STEP[(3, 5)])
    return max(ratios)


def aba_cylinder_shape_error() -> float:
    """The characteristic reference's shape-metric error ESTIMATE at the A|B|A cylinder's working point: its panel step
    up over one minus :func:`largest_measured_ratio`, plus the transport block's step.

    The transport term is the k step (:data:`ABA_CYLINDER_TRANSPORT_STEP`) standing in for a shape step nobody
    measured ``[R]``: 6.9e-11 against a panel term of 2.6e-4, so its exact value cannot move the tolerance.
    """
    p = WORKING_DEGREE["aba_cylinder"]
    return ABA_CYLINDER_SHAPE_STEP[(p, p + 1)] / (1.0 - largest_measured_ratio()) + ABA_CYLINDER_TRANSPORT_STEP


def aba_sphere_shape_error() -> float:
    """The characteristic reference's shape-metric error ESTIMATE at the A|B|A sphere's working point: geometric over
    steps of two rungs, the shape converging in pairs of degrees (:data:`ABA_SPHERE_SHAPE_STEP`)."""
    p = WORKING_DEGREE["aba_sphere"]
    return geometric_error(ABA_SPHERE_SHAPE_STEP[(p, p + 2)], ABA_SPHERE_SHAPE_STEP[(p - 2, p)])


# ── the SN side, relative residuals at the fixtures ───────────────────────
#: ``sphere-sn``: Gauss-Legendre 32, 40 cells (the fixture) against 160 cells (mesh step) and 64 ordinates at 160
#: cells (angular step), k relative and shape metric. ``[M]`` 2026-09-26 (moved here from
#: ``_trajectory_resolvent_ladders`` at step (d): they are the SN side's, and survive the old family).
SPHERE_3REG_SN_STEPS = {
    "k": {"mesh": 5.5e-6, "angle": 9.7e-6},
    "shape": {"mesh": 3.74e-3, "angle": 3.52e-3},
}
#: ``cylinder-sn``: folded 16x32, 40 cells (the fixture) against 32x64 at 40 cells (angular step) and 16x32 at 160
#: cells (mesh step). ``[M]`` 2026-09-26.
CYLINDER_3REG_SN_STEPS = {
    "k": {"mesh": 1.11e-5, "angle": 1.88e-5},
    "shape": {"mesh": 1.59e-3, "angle": 3.44e-3},
}
#: ``cylinder-sn``: the folded 4x8, 40-cell solve's k against its 32x64 limit (1.2310213 against 1.2317490): the
#: unified and standoff rows' SN side. ``[M]`` 2026-09-26.
CYLINDER_3REG_SN_4X8_K_STEP = 5.9e-4


# ── re-measurement ────────────────────────────────────────────────────────


def _timed(label: str, read):
    start = time.perf_counter()
    value = read()
    print(f"{label}: {value!r} ({time.perf_counter() - start:.1f} s)", flush=True)
    return value


def _aba_sphere() -> None:
    """k and the 80 gauged cell averages at degrees 3 to 8, in process; every step between two rungs."""
    from orpheus.derivations.continuous.characteristic import characteristic_reference
    from orpheus.geometry import CoordSystem
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.traced_memo import bypass
    from tests.gates.sn.regression import _generate_snapshots as snapshots
    from tests.gates.sn.verification.analytical import _aba_reference as aba

    observables = aba.shape_observables(snapshots._sphere_3region("2g", 40)["mesh"])
    specification = aba.aba_specification(CoordSystem.SPHERICAL)
    k, shape = {}, {}
    with bypass():
        for p in range(3, 9):
            reference = characteristic_reference(specification, rung(p))
            k[p] = _timed(f"aba sphere p = {p}: k", lambda: reference.read(Eigenvalue()).value)
            shape[p] = [reference.read(ratio).value for _, _, ratio in observables]
    for low in range(3, 8):
        for high in range(low + 1, 9):
            scale = max(abs(v) for v in shape[high])
            step = max(abs(a - b) for a, b in zip(shape[low], shape[high])) / scale
            print(f"aba sphere {low} -> {high}: k step {_relative_step(k, low, high):.4e}; shape step {step:.4e}")


def _aba_cylinder(degrees: tuple[int, ...]) -> None:
    """k and the 80 gauged cell averages on the SN rows' 40-cell mesh at each requested panel degree (default 2, 3, 4,
    :func:`cylinder_panel_rung`), in process, one rung at a time (the p = 3 solve is about 21 min), each rung's shape
    printed whole (``float.hex``) so that the steps can be formed from the log."""
    from orpheus.derivations.continuous.characteristic import characteristic_reference
    from orpheus.geometry import CoordSystem
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.traced_memo import bypass
    from tests.gates.sn.regression import _generate_snapshots as snapshots
    from tests.gates.sn.verification.analytical import _aba_reference as aba

    observables = aba.shape_observables(snapshots._cylinder_3region("2g", 40, "folded_4x8")["mesh"])
    specification = aba.aba_specification(CoordSystem.CYLINDRICAL)
    for p in degrees:
        with bypass():
            reference = characteristic_reference(specification, cylinder_panel_rung(p))
            _timed(f"aba cylinder p = {p}: k", lambda: reference.read(Eigenvalue()).value)
            shape = [reference.read(ratio).value for _, _, ratio in observables]
        print(f"aba cylinder p = {p}: shape {[v.hex() for v in shape]}", flush=True)


def partial_reflector_specification(body: str) -> GeometrySpecification:
    """ERR-094's bodies posed for the reference: mixture A at P0 (``isotropic_mixture``), two groups."""
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.data.materials import Materials
    from orpheus.geometry import StructuredGeometry
    from orpheus.geometry.boundary import AlbedoBoundary, SpecularReturn
    from orpheus.numerics.question import Eigen
    from orpheus.specification.specification import GeometrySpecification
    from tests.gates.sn.verification.analytical._aba_reference import isotropic_mixture

    match body:
        case "slab":
            geometry = StructuredGeometry.slab((0.0, 4.0), (0,), left=AlbedoBoundary(0.3, SpecularReturn("x")),
                                               right=AlbedoBoundary(0.7, SpecularReturn("x")))
        case "sphere":
            geometry = StructuredGeometry.sphere((0.0, 4.0), (0,), outer=AlbedoBoundary(0.7, SpecularReturn("x")))
        case _:
            raise KeyError(body)
    return GeometrySpecification(Materials({0: isotropic_mixture("A")}), geometry,
                                 Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


def edge_specification(body: str) -> GeometrySpecification:
    """The homogeneous edge bodies (the ``*_homogeneous_*_dd_n20`` snapshots' problems): mixture A at P0, R = 2 cm,
    reflective; the sphere at two groups, the cylinder at one."""
    from orpheus.data.cells import CellCoefficient, Channel
    from orpheus.data.materials import Materials
    from orpheus.geometry import BC, CoordSystem, StructuredGeometry
    from orpheus.numerics.question import Eigen
    from orpheus.specification.specification import GeometrySpecification
    from tests.gates.sn.verification.analytical._aba_reference import isotropic_mixture

    coord, groups = {"sphere_2g": (CoordSystem.SPHERICAL, "2g"), "cylinder_1g": (CoordSystem.CYLINDRICAL, "1g")}[body]
    geometry = StructuredGeometry(coord=coord, breakpoints=(0.0, 2.0), mat_ids=(0,), boundaries=(BC.reflective,))
    return GeometrySpecification(Materials({0: isotropic_mixture("A", groups)}), geometry,
                                 Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


def edge_kinf(body: str) -> float:
    """The edge body's medium's k_inf, by the closed form (``kinf_homogeneous``) on the library's arrays."""
    import numpy as np

    from orpheus.derivations.common.eigenvalue import kinf_homogeneous
    from orpheus.derivations.common.xs_library import get_xs

    xs = get_xs("A", {"sphere_2g": "2g", "cylinder_1g": "1g"}[body])
    return float(kinf_homogeneous(np.asarray(xs["sig_t"]), np.asarray(xs["sig_s"]), np.asarray(xs["nu"] * xs["sig_f"]),
                                  np.asarray(xs["chi"])))


def _partial_reflector() -> None:
    from orpheus.derivations.continuous.characteristic import characteristic_reference
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.traced_memo import bypass

    for body in ("slab", "sphere"):
        k = {}
        with bypass():
            for p in range(3, 9):
                reference = characteristic_reference(partial_reflector_specification(body), rung(p))
                k[p] = _timed(f"partial reflector {body} p = {p}: k", lambda: reference.read(Eigenvalue()).value)
        for low in range(3, 8):
            print(f"{body} {low} -> {low + 1}: {_relative_step(k, low, low + 1):.4e}; {low} -> 8: {_relative_step(k, low, 8):.4e}")


def _edge() -> None:
    from orpheus.derivations.continuous.characteristic import characteristic_reference
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.traced_memo import bypass

    import numpy as np

    from tests.gates.sn._test_helpers import SN_TESTS_ROOT

    for snapshot_id, body in _EDGE_SNAPSHOT_BODY.items():
        k = float(np.load(SN_TESTS_ROOT / "regression" / "snapshots" / f"{snapshot_id}.npz")["keff"])
        print(f"snapshot {snapshot_id}: k {k!r}, |k - k_inf| / k_inf = {abs(k - edge_kinf(body)) / edge_kinf(body):.4e}")
    for body in ("sphere_2g", "cylinder_1g"):
        kinf = edge_kinf(body)
        print(f"edge {body}: k_inf {kinf!r}", flush=True)
        with bypass():
            for p in (3, 4):
                reference = characteristic_reference(edge_specification(body), rung(p))
                k = _timed(f"edge {body} p = {p}: k", lambda: reference.read(Eigenvalue()).value)
                print(f"   |k - k_inf| / k_inf = {abs(k - kinf) / kinf:.3e}")


def _sn(geometry: str) -> None:
    """The SN rungs' k, and the shape steps from the fixture in the gauged cell-average metric."""
    from orpheus.numerics.observable import Eigenvalue
    from orpheus.numerics.quadrature import Quadrature
    from tests.gates.sn.regression import _generate_snapshots as generator
    from tests.gates.sn.verification.analytical import _aba_reference as aba
    if geometry == "sphere":
        rungs = {"fixture": (40, 32), "mesh": (160, 32), "angle": (160, 64)}
        build = lambda n, q: {**generator._sphere_3region("2g", n),
                              "quadrature": Quadrature.gauss_legendre(n_ordinates=q), "max_inner": 2000}
    else:
        rungs = {"fixture": (40, (16, 32)), "angle": (40, (32, 64)), "mesh": (160, (16, 32)), "4x8": (40, (4, 8))}
        build = lambda n, q: {**generator._cylinder_3region("2g", n, "folded_4x8"),
                              "quadrature": Quadrature.folded_product(n_mu=q[0], n_phi=q[1]), "max_inner": 2000}
    observables = aba.shape_observables(build(40, rungs["fixture"][1])["mesh"])
    solved = {}
    for name, (n, q) in rungs.items():
        result = generator.run_case(build(n, q))
        k = result.read(Eigenvalue()).value
        solved[name] = (k, {(i, g): result.read(ratio).value for i, g, ratio in observables})
        print(f"{name} {n} {q}: k = {k!r}", flush=True)
    k0, shape0 = solved["fixture"]
    for name in ("mesh", "angle", "4x8"):
        if name in solved:
            k, shape = solved[name]
            scale = max(abs(v) for v in shape.values())
            print(f"{name} step: k {abs(k - k0) / k:.3e}, shape {max(abs(shape0[key] - shape[key]) for key in shape0) / scale:.3e}")


def _records() -> None:
    """The cylinder RECORD's readings (``_aba_reference.CYLINDER_3REG_RECORD``), through the rows' own helpers."""
    from orpheus.numerics.observable import Eigenvalue
    from tests.gates.sn.sweep.curvilinear import test_unified_matvec_cylinder as unified
    from tests.gates.sn.verification.analytical import test_l1_standoff_slab_cylinder as standoff
    from tests.gates.sn.verification.analytical import test_phase_c_crosscheck as pc
    print("phase C:", pc._cylinder_readings(), flush=True)
    print("unified k:", repr(unified._unified_cylinder_solution().read(Eigenvalue()).value), flush=True)
    print("standoff sweep k (nx=40):", repr(standoff._solve_cyl_via_sweep(nx=40).read(Eigenvalue()).value), flush=True)


_LADDERS = {
    "aba-sphere": lambda *_: _aba_sphere(),
    "aba-cylinder": lambda *degrees: _aba_cylinder(tuple(int(d) for d in degrees) or (2, 3, 4)),
    "partial-reflector": lambda *_: _partial_reflector(),
    "edge": lambda *_: _edge(),
    "sphere-sn": lambda *_: _sn("sphere"),
    "cylinder-sn": lambda *_: _sn("cylinder"),
    "records": lambda *_: _records(),
}

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(prog="python -O -m tests.gates.derivations._characteristic_ladders",
                                     description="Re-measure one of the characteristic reference's ladders.")
    parser.add_argument("ladder", choices=sorted(_LADDERS))
    parser.add_argument("degrees", nargs="*", help="aba-cylinder only: the panel degrees (default 2 3 4)")
    arguments = parser.parse_args()
    _LADDERS[arguments.ladder](*arguments.degrees)
