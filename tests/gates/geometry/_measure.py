r"""The one measure of a coordinate system, written in the test from the geometry formulas (layer DAG ruling 31).

The measure of an interval :math:`[a, b]` of a 1-D coordinate system is :math:`m = c\,(T(b) - T(a))`, :math:`T(r) = r^d`:
:math:`d = 1, 2, 3` and :math:`c = 1, \pi, 4\pi/3` for the slab, the cylinder (per unit height) and the sphere. Production
owns the definition (``CoordSystem.measure``, ``orpheus/geometry/coord.py``); this module is the tests' hand-written oracle
of it, and the ONE place the tests write it: ``tests/gates/geometry/test_coord.py`` (the definition pin),
``tests/gates/mesh/test_partition.py`` (the interval rules), ``tests/gates/mesh/test_mesher.py`` (the meshes of named
geometries) and ``tests/gates/mesh/test_mesh1d.py`` (the volume band) import it. Until 2026-10-26 the first two wrote it
twice (``_EXPONENT``, ``_CONSTANT``, ``_measure`` and ``_MEASURE_POWER``, ``_MEASURE_PREFACTOR``,
``_closed_form_region_volume``) and the third wrote the exponent a third time.

It is an ORACLE: it reads no production measure (``CoordSystem.measure``, ``measure_constant``, ``measure_coordinate``,
``compute_volumes_1d``), because a pin routed through production compares the definition with itself
(``retirement-audit`` D.14). ``CoordSystem`` is imported for its members as dictionary keys only. The witness of both
properties is ``tests/gates/geometry/test_coord.py::TestTheMeasureOracleHasOneHome``.

The interval populations and the ``k of N`` reporter live here because the measure laws of ``test_coord.py`` and the rule
laws of ``test_partition.py`` read the same population.
"""
from __future__ import annotations

import numpy as np

from orpheus.geometry import CoordSystem, StructuredGeometry

_SLAB = CoordSystem.CARTESIAN
_CYLINDER = CoordSystem.CYLINDRICAL
_SPHERE = CoordSystem.SPHERICAL

#: The measure's exponent d and constant c, written from geometry here and
#: never read from ``CoordSystem``.
EXPONENT = {_SLAB: 1, _CYLINDER: 2, _SPHERE: 3}
CONSTANT = {_SLAB: 1.0, _CYLINDER: np.pi, _SPHERE: (4.0 / 3.0) * np.pi}


def measure(coord: CoordSystem, edges) -> np.ndarray:
    """The ruled definition, written in the test: ``c * diff(edges**d)`` on arrays."""
    e = np.asarray(edges, dtype=float)
    return CONSTANT[coord] * np.diff(e ** EXPONENT[coord])


def r2_edges(p: int, a: float, b: float, n: int) -> np.ndarray:
    """The ruled R2 body ``T^-1(T(a) + f_j (T(b) - T(a)))``, ends pinned, on arrays.

    The definition pin, not an independent reference (a rule and this copy
    agree if both are transcribed from the ruling); the independent anchors
    are the sum law, the refinement law and the ERR-020 invariant.
    """
    t_a, t_b = np.array([a, b]) ** p
    t = t_a + np.linspace(0.0, 1.0, n + 1) * (t_b - t_a)
    edges = {1: t.copy(), 2: np.sqrt(t), 3: np.cbrt(t)}[p]
    edges[0], edges[-1] = a, b
    return edges


def uniform_geometry(coord: CoordSystem, breakpoints: tuple[float, ...]) -> StructuredGeometry:
    """A geometry with one material per interval and tag 0 on every face (a tag cannot be a deck on a curved face)."""
    return StructuredGeometry.uniform_boundary(
        coord, breakpoints, tuple(range(len(breakpoints) - 1)), 0,
    )


def require(condition: bool, message: str) -> None:
    """A helper's assertion that survives ``python -O``."""
    if not condition:
        raise AssertionError(message)


def report(failures: list[str], population: int, what: str) -> None:
    """Fail with ``k of N`` and the first failures; an empty population is a
    broken harness, never a pass."""
    require(population > 0, f"{what}: the population is empty")
    require(
        not failures,
        f"{what}: {len(failures)} of {population} cases fail; first: " + "; ".join(failures[:5]),
    )


#: The probe's own interval population (S3.6), hollow bodies included.
INTERVALS: tuple[tuple[float, float], ...] = (
    (0.0, 1.0), (0.0, 3.0), (0.0, 2.872), (0.9, 1.1), (1.1, 1.8),
    (0.0, 0.41), (0.01, 2.0), (2.0, 7.0),
)


def seeded_intervals(k: int, seed: int) -> tuple[tuple[float, float], ...]:
    """``k`` intervals, half from the origin, drawn from one fixed seed."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(k):
        a = float(rng.uniform(0.0, 3.0)) if rng.random() < 0.5 else 0.0
        out.append((a, a + float(rng.uniform(0.01, 5.0))))
    return tuple(out)


#: ``(0.3, 0.7)``: on a sphere ``c * (b**3 - a**3)`` and ``c * b**3 - c * a**3``
#: differ there, so a re-association of the measure is visible (spec §1.3a,
#: premise 6).
LAW_INTERVALS = INTERVALS + ((0.3, 0.7),) + seeded_intervals(8, seed=405)
