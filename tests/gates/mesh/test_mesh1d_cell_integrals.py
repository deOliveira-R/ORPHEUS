r"""The cell co-vector of a mesh-free weight on a 1-D mesh (#405 P2 step 7b.2.1, gates R7b2.8).

``Mesh1D.cell_integrals(weight) -> (G, N)``: entry ``[g, i]`` is
:math:`\int_{\text{cell } i} w_g\,dV` in the coordinate system's measure
(:math:`dV = c\,T'(r)\,dr`: ``dr`` on a slab, :math:`2\pi r\,dr` on a
cylinder per unit height, :math:`4\pi r^2\,dr` on a sphere), the layout of
``ScalarFlux.values`` on a 1-D mesh. It is the mesh's realisation of a
:class:`~orpheus.numerics.observable.FluxIntegral`'s weight, so a cell-average
answer reads the flux integral exactly as :math:`\sum_{g,i}\phi_{g,i}
\,\texttt{cell\_integrals}(w)_{g,i}` (spec §1.7b.2, "7b.2.1").

* A ``RegionwiseConstant`` ``(R, G)`` reads through the region labels
  (7b.2.0): ``per_cell(values).T * volumes``; ``R`` must be the mesh's region
  count.
* A ``Symbolic`` weight is integrated EXACTLY by SymPy in ``Symbolic.r``
  over each cell with the measure density, then rounded to a float; a
  ``Piecewise`` whose step lies inside a cell integrates exactly; a weight
  that depends on the direction is refused (``Symbolic.depends_on``, the one
  definition).

The references below are closed forms of the measure, written independently
of the production formula (never ``per_cell(...) * volumes`` again). Gate ids
``R7b2.8.<n>``. First red on ``fa38de31``: ``AttributeError``, ``Mesh1D`` has
no ``cell_integrals``. Every test is ``foundation``, claim kind THEOREM.
"""
from __future__ import annotations

import math

import numpy as np
import pytest
import sympy

from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_mesh1d_cell_integrals.py"
_LABELS = "tests/gates/mesh/test_mesh1d_regions.py::test_r7b2_0_4_aba_regions_0_and_2_are_distinct"
_SLAB, _CYLINDER, _SPHERE = CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL
_COORDS = (_SLAB, _CYLINDER, _SPHERE)

#: The measure of [a, b] in closed form, written from the geometry (not read off the mesh).
_SHELL = {
    _SLAB: lambda a, b: b - a,
    _CYLINDER: lambda a, b: math.pi * (b * b - a * a),
    _SPHERE: lambda a, b: 4.0 / 3.0 * math.pi * (b**3 - a**3),
}
#: The integral of r over [a, b] in the same measure, closed form (catches an exponent error the constant hides).
_FIRST_MOMENT = {
    _SLAB: lambda a, b: (b * b - a * a) / 2.0,
    _CYLINDER: lambda a, b: 2.0 * math.pi * (b**3 - a**3) / 3.0,
    _SPHERE: lambda a, b: math.pi * (b**4 - a**4),
}


def _aba(coord: CoordSystem, cells=(2, 4, 2)) -> Mesh1D:
    """A|B|A at 0.5, 1.5, 2.0 cm, materials (0, 1, 0); cells of width 0.25 cm."""
    g = StructuredGeometry.from_thicknesses(
        coord=coord, thicknesses=(0.5, 1.0, 0.5), mat_ids=(0, 1, 0),
        boundaries=(BC.vacuum,) if coord is not _SLAB else (BC.reflective, BC.vacuum),
    )
    return Mesher(g).partition(tuple(CellsByCount.uniform_width(n) for n in cells)).mesh


def _close(actual: float, expected: float, ulps: int = 8) -> bool:
    """Within ``ulps`` ulp of the larger magnitude: an exact integral rounded once, against a closed form in floats."""
    return abs(actual - expected) <= ulps * math.ulp(max(abs(actual), abs(expected), 1e-300))


@pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
def test_r7b2_8_1_the_layout_and_the_indicator_of_each_cell(coord) -> None:
    """``(G, N)``; the weight 1 in group 0 and 0 in group 1 gives each cell's measure in row 0 and zeros in row 1."""
    mesh = _aba(coord)
    out = mesh.cell_integrals(Symbolic.of(1, 0))
    assert out.shape == (2, mesh.N)
    edges = np.asarray(mesh.edges)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        assert _close(float(out[0, i]), _SHELL[coord](a, b)), (i, out[0, i], _SHELL[coord](a, b))
    np.testing.assert_array_equal(out[1], np.zeros(mesh.N))


@pytest.mark.rests_on(f"{_HERE}::test_r7b2_8_1_the_layout_and_the_indicator_of_each_cell")
@pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
def test_r7b2_8_2_the_sum_over_cells_is_the_geometry_volume(coord) -> None:
    """Σ_i ∫ 1 dV is the body's measure, 4π R³/3, π R² or R (closed form), to N·8 ulp."""
    mesh = _aba(coord)
    total = float(np.sum(mesh.cell_integrals(Symbolic.of(1, 1))[0]))
    assert _close(total, _SHELL[coord](0.0, 2.0), ulps=8 * mesh.N)


@pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
def test_r7b2_8_3_the_measure_density(coord) -> None:
    """The weight r reads the first moment in each coordinate system's measure (b⁴ − a⁴ on a sphere): an exponent
    or a constant taken from another coordinate system reddens this row and not only the volume rows."""
    mesh = _aba(coord)
    out = mesh.cell_integrals(Symbolic.of(Symbolic.r, Symbolic.r))
    edges = np.asarray(mesh.edges)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        assert _close(float(out[1, i]), _FIRST_MOMENT[coord](a, b)), (i, out[1, i])


def test_r7b2_8_3_sphere_and_cylinder_measures_differ() -> None:
    """Activation of the coordinate system: one weight on equal edges reads differently on a sphere and a cylinder."""
    sphere, cylinder = _aba(_SPHERE).cell_integrals(Symbolic.of(1, 1)), _aba(_CYLINDER).cell_integrals(Symbolic.of(1, 1))
    assert not np.any(sphere == cylinder)


@pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
def test_r7b2_8_4_a_step_inside_a_cell_integrates_exactly(coord) -> None:
    """The weight 1 on r < 0.6, a step inside the cell [0.5, 0.75]: that cell reads the partial measure of
    [0.5, 0.6], the cells below read their full measure, the cells above read 0. A midpoint or nodal rule
    would read the cell's whole measure or none of it."""
    mesh = _aba(coord)
    step = sympy.Piecewise((1, Symbolic.r < sympy.Rational(3, 5)), (0, True))
    out = mesh.cell_integrals(Symbolic.of(step, 0))[0]
    edges = np.asarray(mesh.edges)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        expected = _SHELL[coord](a, b) if b <= 0.6 else (_SHELL[coord](a, 0.6) if a < 0.6 else 0.0)
        assert _close(float(out[i]), expected), (i, out[i], expected)
    assert 0.0 < out[2] < _SHELL[coord](0.5, 0.75)  # activation: the straddling cell is the third


@pytest.mark.rests_on(_LABELS)
@pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
def test_r7b2_8_5_a_region_table_reads_through_the_labels(coord) -> None:
    """A ``RegionwiseConstant`` with rows 0 and 2 different reads them separately on the two fuel regions of
    A|B|A (one material, two regions: the X1 witness 7b.2.0 enabled). The reference finds each cell's region
    by containment of its centre in the breakpoints, independently of the labels."""
    mesh = _aba(coord)
    table = np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])  # (regions, groups)
    out = mesh.cell_integrals(RegionwiseConstant(table))
    centres = 0.5 * (np.asarray(mesh.edges)[:-1] + np.asarray(mesh.edges)[1:])
    region = np.searchsorted([0.5, 1.5, 2.0], centres)
    edges = np.asarray(mesh.edges)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        for g in range(2):
            assert _close(float(out[g, i]), table[region[i], g] * _SHELL[coord](a, b)), (g, i)
    fuel = mesh.mat_ids == 0
    assert set(np.round(out[0, fuel] / np.asarray(mesh.volumes)[fuel], 12).tolist()) == {1.0, 3.0}


def test_r7b2_8_6_refusals() -> None:
    """A region table of the wrong region count; a weight that depends on the direction (refused by the one
    definition, ``Symbolic.depends_on``)."""
    mesh = _aba(_SPHERE)
    with pytest.raises(ValueError, match="region"):
        mesh.cell_integrals(RegionwiseConstant(np.ones((2, 2))))
    for direction in (Symbolic.mu, Symbolic.phi):
        with pytest.raises(ValueError, match="depends on"):
            mesh.cell_integrals(Symbolic.of(direction, 0))


def test_r7b2_8_6_the_direction_rule_is_the_one_definition(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE: the refusal consults ``Symbolic.depends_on`` (rebinding it to say "independent" admits the weight)."""
    mesh = _aba(_SLAB, (1, 1, 1))
    calls: list[tuple] = []

    def never(self, *coordinates):
        calls.append(coordinates)
        return False

    monkeypatch.setattr(Symbolic, "depends_on", never)
    try:
        mesh.cell_integrals(Symbolic.of(Symbolic.mu + 1, 0))
    except ValueError as error:
        assert "depends on the direction" not in str(error), "a second direction rule refused it"
    except Exception:  # noqa: BLE001 - the integrator may fail on mu; the route is what is asserted
        pass
    assert calls, "Symbolic.depends_on was never consulted"


# ─────────────────────────────────────────────────────────────────────
# R7b2.8.7 — the scope edge of the exact integration (the review round of 7b.2.1)
# ─────────────────────────────────────────────────────────────────────


def _quad_reference(coord: CoordSystem, f, a: float, b: float, points) -> float:
    """∫_a^b f(r) dV by adaptive quadrature split at the given points: a route independent of SymPy."""
    from scipy.integrate import quad

    density = {_SLAB: lambda r: 1.0, _CYLINDER: lambda r: 2.0 * math.pi * r, _SPHERE: lambda r: 4.0 * math.pi * r * r}[coord]
    inside = [p for p in points if a < p < b]
    return quad(lambda r: f(r) * density(r), a, b, points=inside or None, epsabs=0.0, epsrel=1e-13, limit=200)[0]


def test_r7b2_8_7_a_step_on_a_non_polynomial_argument_is_refused() -> None:
    """SCOPE-BOUNDARY: ``Piecewise((1, sin(3r) > 0), (0, True))`` is refused (``[M]`` qa: SymPy 1.14 integrates
    it over [2, 3] as 0); the refusal names the boundary."""
    mesh = _aba(_SPHERE)
    step = sympy.Piecewise((1, sympy.sin(3 * Symbolic.r) > 0), (0, True))
    with pytest.raises(ValueError, match="is integrated only where its argument is polynomial in r"):
        mesh.cell_integrals(Symbolic.of(step, 0))


@pytest.mark.parametrize("coord", _COORDS, ids=lambda c: c.name.lower())
def test_r7b2_8_7_polynomial_arguments_are_admitted_and_exact(coord) -> None:
    """Inside the boundary: a step on ``r**2 < 2`` (its edge √2 inside the cell [1.25, 1.5]) and ``Abs(r - 1.1)``
    (its kink inside [1.0, 1.25]) integrate exactly, against adaptive quadrature split at the kink (an independent
    route), to 1e-12 relative."""
    mesh = _aba(coord)
    r = Symbolic.r
    cases = (
        (sympy.Piecewise((1, r**2 < 2), (0, True)), lambda x: 1.0 if x * x < 2 else 0.0, (math.sqrt(2.0),)),
        (sympy.Abs(r - sympy.Rational(11, 10)), lambda x: abs(x - 1.1), (1.1,)),
    )
    edges = np.asarray(mesh.edges)
    for expression, f, kinks in cases:
        out = mesh.cell_integrals(Symbolic.of(expression, 0))[0]
        for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
            expected = _quad_reference(coord, f, float(a), float(b), kinks)
            assert abs(float(out[i]) - expected) <= 1e-12 * max(abs(expected), _SHELL[coord](a, b)), (expression, i)


def test_r7b2_8_7_an_unintegrable_weight_gets_the_keyed_refusal() -> None:
    """``floor(2r)`` (a bare ``TypeError`` inside SymPy before the review round) is a keyed ``ValueError``."""
    mesh = _aba(_SLAB, (1, 1, 1))
    with pytest.raises(ValueError, match="cannot be integrated over the cells"):
        mesh.cell_integrals(Symbolic.of(sympy.floor(2 * Symbolic.r), 0))


# ─────────────────────────────────────────────────────────────────────
# R7b2.8.8 — a hollow body (r_0 > 0)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("coord", [_CYLINDER, _SPHERE], ids=lambda c: c.name.lower())
def test_r7b2_8_8_a_hollow_body(coord) -> None:
    """A shell starting at r_0 = 0.5: the indicator and the first moment of each cell against the closed forms,
    and the sum against the shell's measure (an antiderivative evaluated from 0 rather than r_0 reddens)."""
    g = StructuredGeometry(coord=coord, breakpoints=(0.5, 1.0, 2.0), mat_ids=(0, 1), boundaries=(BC.reflective, BC.vacuum))
    mesh = Mesher(g).partition(CellsByCount.uniform_width(2)).mesh
    ones = mesh.cell_integrals(Symbolic.of(1))[0]
    moments = mesh.cell_integrals(Symbolic.of(Symbolic.r))[0]
    edges = np.asarray(mesh.edges)
    for i, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        assert _close(float(ones[i]), _SHELL[coord](a, b)) and _close(float(moments[i]), _FIRST_MOMENT[coord](a, b)), i
    assert _close(float(np.sum(ones)), _SHELL[coord](0.5, 2.0), ulps=8 * mesh.N)
