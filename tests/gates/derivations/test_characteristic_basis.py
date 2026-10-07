"""Gates for the characteristic reference's panel basis (:class:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis`).

P1 step (b), second rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b),
second rung: API sketch", item 1, and its "Ruled 2026-10-06" line).
Re-specification: ``scratch/characteristic_architecture/p1_step_b2/spec.md``,
rows P1-P9 (the spec's C12 and C6b, and the mass rows), plus the re-keyed
walls (W1, W2).

The ladder, bottom up (``rests_on`` on each row):

1. the kernel's partition (``tests/gates/geometry/test_chord.py``) and the
   reference kernel's Gauss rules (``composite_gauss_legendre``);
2. the panel ends: every breakpoint is one, each panel lies in one region
   [P1]; the grading law toward walls and interfaces and not toward a singular
   stratum [P2];
3. the nodes against Gauss-Legendre roots in mpmath [P3]; the functions are
   cardinal and a partition of unity [P4]; the interpolant reproduces every
   per-region polynomial of degree p and no polynomial of degree p + 1 [P5,
   the spec's C12];
4. the metric: the volume density against the hand-written densities [P6];
   the mass matrix against mpmath volume integrals of Lagrange products [P7];
   each panel's mass is its chart measure [P8];
5. the basis reads only the body's partition and the resolution [P9, C6b];
6. the walls re-keyed onto the panel partition [W1, W2].

Every reference is written in this file or in ``_characteristic_mp.py``
(mpmath; nothing there imports the code under test). The panel ends enter the
mass reference as DATA (the panels the basis is defined on), gated on their own
by P1 and P2.

Every row is ``foundation``: the theory page carries no label for the basis
(``characteristic-traversal-integrals`` names tau_k and B_k, which the transport
file's rows verify).
"""
from __future__ import annotations

import inspect

import mpmath as mp
import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic import Walls
from orpheus.derivations.continuous.characteristic.basis import PanelBasis
from orpheus.geometry.boundary import AlbedoBoundary, PeriodicBoundary, SpecularReturn, VacuumInflow
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.structured_geometry import StructuredGeometry

from . import _characteristic_mp as R

_HERE = "tests/gates/derivations/test_characteristic_basis.py::"
_CHORD = "tests/gates/geometry/test_chord.py::"
_WALLS = "tests/gates/derivations/test_characteristic_walls.py::"

_EPS = float(np.finfo(float).eps)

# Spec §4 fixtures (non-dyadic radii): MR3 solid, MR3H hollow, SLB3 slab; plus
# the afterthoughts: a 1e3 width ratio between neighbouring regions and a
# one-region body.
_MR3 = (0.0, 0.5, 1.5, 2.0)
_MR3H = (0.4, 0.5, 1.5, 2.0)
_SLB3 = (0.0, 0.4, 1.5, 2.3)
_THIN_THICK = (0.0, 1.3e-3, 1.3, 1.7)
_ONE = (0.0, 1.3)


def _spec(a: float):
    return AlbedoBoundary(a, SpecularReturn(axis="x"))


def _bodies() -> dict[str, StructuredGeometry]:
    return {
        "sphere_solid": StructuredGeometry.sphere(_MR3, (0, 1, 2), outer=_spec(0.6)),
        "cylinder_solid": StructuredGeometry.cylinder(_MR3, (0, 1, 2), outer=_spec(0.6)),
        "sphere_hollow": StructuredGeometry.sphere(_MR3H, (0, 1, 2), inner=_spec(0.3), outer=_spec(0.8)),
        "cylinder_hollow": StructuredGeometry.cylinder(_MR3H, (0, 1, 2), inner=_spec(0.3), outer=_spec(0.8)),
        "slab": StructuredGeometry.slab(_SLB3, (0, 1, 2), left=_spec(0.3), right=_spec(0.8)),
        "slab_thin_thick": StructuredGeometry.slab(_THIN_THICK, (0, 1, 2), left=VacuumInflow(), right=VacuumInflow()),
        "sphere_one_region": StructuredGeometry.sphere(_ONE, (0,), outer=VacuumInflow()),
        "cylinder_thin_thick": StructuredGeometry.cylinder(_THIN_THICK, (0, 1, 2), outer=VacuumInflow()),
    }


_BODIES = _bodies()
_CHART = {"sphere": "sphere", "cylinder": "cylinder", "slab": "slab"}

#: (degree, layers, ratio): the rung's working point, no grading, a deep grading, piecewise constants.
_RESOLUTIONS = [(3, 2, 0.5), (1, 0, 0.5), (5, 4, 0.3), (0, 1, 0.5)]


def _chart_name(name: str) -> str:
    return name.split("_")[0]


def _basis(geometry: StructuredGeometry, degree: int, layers: int, ratio: float) -> PanelBasis:
    return PanelBasis.of(ConcentricPartition.of(geometry), degree, layers, ratio)


def _cases():
    return [(f"{name}-p{p}L{lay}q{q}", name, p, lay, q) for name in _BODIES for (p, lay, q) in _RESOLUTIONS]


_CASES = _cases()
_IDS = [c[0] for c in _CASES]
_ARGS = [c[1:] for c in _CASES]


def _conditioning(a: float, b: float) -> float:
    """kappa = max(|a|, |b|) / (b - a): the digits a panel's local coordinate loses to its global position.

    The panel's map to [-1, 1] reads ``(2c - (a + b)) / (b - a)``, whose
    rounding is eps |c| / (b - a) relative to the reference interval: a thin
    graded panel far from the origin (width 5e-6 at r = 2) loses five digits
    there whatever the code does, so every per-panel tolerance below carries
    this factor (measured: the worst residual over (1 + kappa) is O(10) ulp).
    """
    return max(abs(a), abs(b)) / (b - a)


def _solid(name: str) -> bool:
    return _BODIES[name].breakpoints[0] == 0.0 and not name.startswith("slab")


# ── P1: the panel ends ────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree", "layers", "ratio"), _ARGS, ids=_IDS)
@pytest.mark.rests_on(_CHORD + "test_every_slot_length_matches_the_closed_form")
def test_every_breakpoint_is_a_panel_end_and_each_panel_lies_in_one_region(name, degree, layers, ratio) -> None:
    """[P1, C12's partition half] The panel ends refine the body's breakpoints with the same ends.

    Every breakpoint (wall and interface) is a panel end, bit for bit; the ends
    are the body's; each panel [a, b] lies inside the region ``region_of_panel``
    names (r_k <= a < b <= r_{k+1}, checked against the breakpoints, not against
    the code's own region table); every region has at least one panel; the
    partition is posed on the body's chart. First red: the partition built over
    [r_0, r_n] without the interior breakpoints (a panel straddling an
    interface).
    """
    geometry = _BODIES[name]
    basis = _basis(geometry, degree, layers, ratio)
    r = geometry.breakpoints
    ends = basis.partition.breakpoints
    missing = [v for v in r if v not in ends]
    assert not missing, f"breakpoints {missing} are not panel ends: {ends}"
    assert ends[0] == r[0] and ends[-1] == r[-1]
    assert basis.partition.chart == Chart(geometry.coord)
    owner = basis.region_of_panel
    assert owner.shape == (len(ends) - 1,)
    for panel, (a, b) in enumerate(zip(ends[:-1], ends[1:])):
        k = int(owner[panel])
        assert r[k] <= a < b <= r[k + 1], f"panel {panel} [{a}, {b}] is not inside region {k} [{r[k]}, {r[k + 1]}]"
    assert sorted(set(owner.tolist())) == list(range(len(r) - 1))


# ── P2: the grading law ───────────────────────────────────────────────────


def _hand_ends(r, solid: bool, layers: int, ratio: float) -> list[R.Mpf]:
    """The panel ends the sketch's grading law places, in mpmath.

    In region [r_k, r_{k+1}], each end that is a wall or an interface (every
    breakpoint except a solid body's r_0 = 0) gets interior ends at distance
    w ratio^j, j = 1..layers, with w the half width when both ends are graded
    and the whole width when one is.
    """
    ends = [mp.mpf(r[0])]
    for k in range(len(r) - 1):
        a, b = mp.mpf(r[k]), mp.mpf(r[k + 1])
        lo = not (solid and k == 0)
        hi = True
        w = (b - a) / 2 if lo and hi else b - a
        inner = []
        if lo:
            inner += [a + w * mp.mpf(ratio) ** j for j in range(layers, 0, -1)]
        if hi:
            inner += [b - w * mp.mpf(ratio) ** j for j in range(1, layers + 1)]
        ends += sorted(inner) + [b]
    return ends


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree", "layers", "ratio"), _ARGS, ids=_IDS)
@pytest.mark.rests_on(_HERE + "test_every_breakpoint_is_a_panel_end_and_each_panel_lies_in_one_region")
def test_panels_grade_toward_walls_and_interfaces_and_not_toward_a_singular_stratum(name, degree, layers, ratio) -> None:
    """[P2] The panel ends are the grading law's, written by hand in mpmath; 4 ulp of the body's size.

    The law (the sketch's item 1, as built): ``layers`` ends at distance
    w ratio^j from each graded end. A solid body's centre or axis (r_0 = 0, a
    singular stratum) is NOT graded; a hollow body's inner wall and every slab
    face ARE. A second, law-free leg: in a solid body's region 0 the panel
    touching the centre is wider than the panel touching r_1 (graded), and in a
    hollow body's region 0 the panel touching the inner wall is the narrowest
    of the region (the control: the inner wall is graded).

    First reds: (a) the centre graded like a wall: the hand set and the
    law-free leg red on every solid body with layers > 0; (b) the inner wall of
    a hollow body left ungraded: the control leg reds; (c) the depths measured
    from the wrong end (w (1 - ratio^j)).
    """
    geometry = _BODIES[name]
    basis = _basis(geometry, degree, layers, ratio)
    solid = _solid(name)
    want = [float(v) for v in _hand_ends(geometry.breakpoints, solid, layers, ratio)]
    got = np.asarray(basis.partition.breakpoints)
    assert got.shape == (len(want),), f"{len(got)} ends, the law places {len(want)}"
    np.testing.assert_allclose(got, want, rtol=0.0, atol=4 * _EPS * abs(geometry.breakpoints[-1]))
    if layers == 0:
        return
    widths = np.diff(got)
    region0 = widths[basis.region_of_panel == 0]
    if solid and (layers > 1 or ratio < 0.5):     # L = 1 at ratio 1/2 cuts region 0 into two equal halves
        assert region0[0] > region0[-1], f"the centre is graded: widths {region0}"
    elif not solid:
        assert region0[0] <= region0.min() * (1.0 + 1e-9), f"the inner wall is not graded: widths {region0}"


# ── P3: the nodes ─────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree", "layers", "ratio"), _ARGS, ids=_IDS)
@pytest.mark.rests_on(_HERE + "test_every_breakpoint_is_a_panel_end_and_each_panel_lies_in_one_region")
def test_the_nodes_are_each_panels_gauss_legendre_points(name, degree, layers, ratio) -> None:
    """[P3] Node m of panel P is the m-th root of P_{p+1} mapped onto the panel (mpmath roots), index P(p+1) + m.

    Also: ``panel`` and ``region`` of each node agree with the panel ends and
    the breakpoints (the node lies strictly inside its panel, never on an end,
    so no node sits on an interface), ``N = P (p + 1)``. 4 ulp of the panel's
    outer end. First reds: nodes of degree p (one node short); a panel's nodes
    mapped onto its neighbour; the panel-major order broken (node-major).
    """
    geometry = _BODIES[name]
    basis = _basis(geometry, degree, layers, ratio)
    ends = basis.partition.breakpoints
    n_panels = len(ends) - 1
    nodes = np.asarray(basis.nodes)
    assert nodes.shape == (n_panels * (degree + 1),)
    want = np.array([float(x) for a, b in zip(ends[:-1], ends[1:]) for x in R.gl_nodes(a, b, degree + 1)])
    scale = np.repeat([max(abs(a), abs(b)) for a, b in zip(ends[:-1], ends[1:])], degree + 1)
    assert np.all(np.abs(nodes - want) <= 4 * _EPS * scale), np.max(np.abs(nodes - want) / scale) / _EPS
    panel = np.asarray(basis.panel)
    np.testing.assert_array_equal(panel, np.repeat(np.arange(n_panels), degree + 1))
    lo, hi = np.asarray(ends)[panel], np.asarray(ends)[panel + 1]
    assert np.all((lo < nodes) & (nodes < hi)), "a node lies on a panel end"
    r = np.asarray(geometry.breakpoints)
    region = np.asarray(basis.region)
    assert np.all((r[region] < nodes) & (nodes < r[region + 1]))


# ── P4: cardinal functions, a partition of unity ──────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree", "layers", "ratio"), _ARGS, ids=_IDS)
@pytest.mark.rests_on(_HERE + "test_the_nodes_are_each_panels_gauss_legendre_points")
def test_the_functions_are_cardinal_at_the_nodes_and_sum_to_one(name, degree, layers, ratio) -> None:
    """[P4] values(node_m of P, P) is the m-th unit vector (16 (1 + kappa) ulp, ``_conditioning``) and the p + 1 functions sum to 1 everywhere on P.

    First reds: the reference nodes of degree p - 1 used for the functions; the
    affine map to [-1, 1] reversed (a palindromic node set makes the CARDINAL
    leg blind to a reversal: declared; the reproduction row P5 sees it through
    the odd powers).
    """
    basis = _basis(_BODIES[name], degree, layers, ratio)
    nodes = np.asarray(basis.nodes)
    panel = np.asarray(basis.panel)
    values = basis.values(nodes, panel)
    assert values.shape == (len(nodes), degree + 1)
    ends = np.asarray(basis.partition.breakpoints)
    kappa = np.array([_conditioning(a, b) for a, b in zip(ends[:-1], ends[1:])])[panel]
    eye = np.tile(np.eye(degree + 1), (len(ends) - 1, 1))
    assert np.all(np.abs(values - eye) <= 16 * _EPS * (1.0 + kappa)[:, None]), \
        float(np.max(np.abs(values - eye) / (1.0 + kappa)[:, None])) / _EPS
    rng = np.random.default_rng(17)
    for p in range(len(ends) - 1):
        c = ends[p] + (ends[p + 1] - ends[p]) * np.concatenate([[0.0, 1.0], rng.random(9)])
        np.testing.assert_allclose(basis.values(c, np.full(c.shape, p)).sum(axis=-1), 1.0, rtol=0.0,
                                   atol=16 * _EPS * (degree + 1))


# ── P5: the exactness family (C12) ────────────────────────────────────────


def _random_polynomial(n_regions: int, degree: int, seed: int, *, even_first: bool = False) -> list[list[float]]:
    """A random polynomial of ``degree`` per region; with ``even_first`` region 0's odd coefficients are 0.

    Region 0 of a solid body holds the even panel (Lagrange in c^2) beside
    ordinary panels, so the region's common exactness family is the EVEN
    polynomials of degree <= p (the spec's RB1).
    """
    rng = np.random.default_rng(seed)
    out = [list(rng.uniform(-1.0, 1.0, degree + 1) + (2.0 * k + 0.5) * (np.arange(degree + 1) == 0))
           for k in range(n_regions)]
    if even_first:
        out[0] = [a if m % 2 == 0 else 0.0 for m, a in enumerate(out[0])]
    return out


def _interpolant_error(basis: PanelBasis, coefficients, n_regions: int, *, conditioned: bool = True) -> float:
    f = R.region_polynomial(coefficients)
    x = R.sample(f, np.asarray(basis.region), np.asarray(basis.nodes))
    rng = np.random.default_rng(5)
    ends = np.asarray(basis.partition.breakpoints)
    worst, scale = 0.0, 0.0
    p1 = basis.degree + 1
    for p in range(len(ends) - 1):
        c = ends[p] + (ends[p + 1] - ends[p]) * np.concatenate([[0.0, 1.0, 1e-9, 1.0 - 1e-9], rng.random(7)])
        got = basis.values(c, np.full(c.shape, p)) @ x[p * p1:(p + 1) * p1]
        k = int(basis.region_of_panel[p])
        want = np.array([float(f(k, R.mpf(v))) for v in c])
        kappa = 1.0 + _conditioning(ends[p], ends[p + 1]) if conditioned else 1.0
        worst = max(worst, float(np.max(np.abs(got - want))) / kappa)
        scale = max(scale, float(np.max(np.abs(want))))
    return worst / scale


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree", "layers", "ratio"), _ARGS, ids=_IDS)
@pytest.mark.rests_on(_HERE + "test_the_functions_are_cardinal_at_the_nodes_and_sum_to_one")
def test_the_interpolant_reproduces_every_per_region_polynomial_of_degree_p(name, degree, layers, ratio) -> None:
    """[P5, C12] Sampling a per-region polynomial of degree p at the nodes and summing the functions returns it.

    The polynomial jumps at every interface (distinct random coefficients per
    region); it is read on every panel at both ends, a hair inside them and at
    seven random points, against mpmath. 64 ulp of the polynomial's maximum
    (the Lebesgue constant of Gauss-Legendre nodes times the evaluation's
    rounding). The loading control: a polynomial of degree p + 1 is NOT
    reproduced (above 1e-8 relative, raw: 1e-7 at the deepest grading, three
    orders above the honest band there), so the row reads the degree. First reds:
    a panel straddling an interface (the jump smeared); the functions of degree
    p - 1; the reference interval's map reversed.

    Re-posed 2026-10-06 (rung 3, the spec's RB1): on a solid body region 0's
    polynomial is EVEN (its odd coefficients zeroed), the family the even panel
    at the stratum and the region's other panels share; the even panel's own
    family (degree 2p in c) is ``test_the_even_panel_spans_the_polynomials_in_c_squared``.
    """
    geometry = _BODIES[name]
    basis = _basis(geometry, degree, layers, ratio)
    n = len(geometry.breakpoints) - 1
    assert _interpolant_error(basis, _random_polynomial(n, degree, 11, even_first=_solid(name)), n) <= 64 * _EPS
    loaded = _interpolant_error(basis, _random_polynomial(n, degree + 1, 13), n, conditioned=False)
    assert loaded > 1e-8, f"a degree-{degree + 1} polynomial is reproduced to {loaded:.1e}: the row cannot see the degree"


# ── P6-P8: the metric ─────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.parametrize("name", list(_BODIES))
def test_the_volume_density_is_the_charts_by_hand(name) -> None:
    """[P6] The basis's chart's ``measure_density(c)`` = 1, 2 pi c, 4 pi c^2 (slab, cylinder per unit height, sphere), written here; 4 ulp.

    Re-pointed 2026-10-06 (rung 3, the spec's RB3): ``PanelBasis.volume_density``
    retired onto the kernel's ``Chart.measure_density`` (the derivative of the
    one measure, gated in ``tests/gates/geometry/test_measure_density.py``);
    the mass reads it through the basis's own chart, which this row reads, and
    the retirement's route (no second density on the basis) is
    ``test_characteristic_assembly.py::test_the_mass_and_the_wall_area_read_the_one_density``.
    First reds: the density measure_constant c^(d-1) (the factor d dropped);
    measure_constant d c^d (one power too many).
    """
    basis = _basis(_BODIES[name], 2, 1, 0.5)
    c = np.array([0.0, 0.37, 1.1, 1.9])
    want = np.array([float(R.density(_chart_name(name), R.mpf(v))) for v in c])
    np.testing.assert_allclose(basis.regions.chart.measure_density(c), want, rtol=4 * _EPS, atol=0.0)


def _mass_reference(chart: str, a: float, b: float, degree: int) -> np.ndarray:
    """The panel's mass block in mpmath; on a radial chart's panel starting at 0 (the stratum) the Lagrange functions are in c^2."""
    nodes = R.gl_nodes(a, b, degree + 1)
    lagrange = R.lagrange_even if (chart != "slab" and a == 0.0) else R.lagrange
    with mp.workdps(R.DPS):
        return np.array([[float(R.quad(lambda c: lagrange(nodes, i, c) * lagrange(nodes, j, c)
                                       * R.density(chart, c), [R.mpf(a), R.mpf(b)]))
                          for j in range(degree + 1)] for i in range(degree + 1)])


_MASS_CASES = [c for c in _CASES if c[2] in (3, 0) or c[1].startswith("sphere_hollow")]


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree", "layers", "ratio"), [c[1:] for c in _MASS_CASES],
                         ids=[c[0] for c in _MASS_CASES])
@pytest.mark.rests_on(_HERE + "test_the_nodes_are_each_panels_gauss_legendre_points",
                      _HERE + "test_the_volume_density_is_the_charts_by_hand")
def test_the_mass_matrix_is_the_volume_integral_of_each_product(name, degree, layers, ratio) -> None:
    """[P7] W_ij = int u_i u_j dV: each panel's block against mpmath, every other entry exactly 0.

    The reference builds its own Lagrange functions on mpmath Gauss-Legendre
    roots over the basis's panel ends and integrates with the density written
    in this file (sphere 4 pi c^2, cylinder 2 pi c, slab 1). 64 (1 + kappa) ulp
    of the block's largest entry (``_conditioning``). First reds: the density without the factor d; the
    mass rule short of exact (re-posed 2026-10-06, rung 3, the spec's RB2: on
    the even panel at a solid body's stratum the integrand has degree 4p + d - 1,
    so the rule takes 2p + 2 points; the rung-2 rule of p + 2 points reds there);
    a block placed at the wrong panel's indices; the weights of one panel scaled
    by its neighbour's width. On the stratum panel the reference's Lagrange
    functions are in c^2 (``R.lagrange_even``, written independently of the basis).
    """
    chart = _chart_name(name)
    basis = _basis(_BODIES[name], degree, layers, ratio)
    W = np.asarray(basis.mass)
    p1 = degree + 1
    ends = basis.partition.breakpoints
    n_panels = len(ends) - 1
    assert W.shape == (n_panels * p1, n_panels * p1)
    block = np.kron(np.eye(n_panels, dtype=bool), np.ones((p1, p1), dtype=bool))
    np.testing.assert_array_equal(W[~block], 0.0)
    for p, (a, b) in enumerate(zip(ends[:-1], ends[1:])):
        got = W[p * p1:(p + 1) * p1, p * p1:(p + 1) * p1]
        want = _mass_reference(chart, a, b, degree)
        scale = float(np.max(np.abs(want))) * (1.0 + _conditioning(a, b))
        assert np.max(np.abs(got - want)) <= 64 * _EPS * scale, (p, np.max(np.abs(got - want)) / scale / _EPS)


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree", "layers", "ratio"), _ARGS, ids=_IDS)
@pytest.mark.rests_on(_HERE + "test_the_functions_are_cardinal_at_the_nodes_and_sum_to_one")
def test_each_panels_mass_is_its_chart_measure(name, degree, layers, ratio) -> None:
    """[P8] sum_ij W_ij over a panel's block is Chart.measure of the panel (the one definition), and the total is the body's.

    The functions sum to 1, so the block's sum is the panel's volume: the basis's
    derived density integrates to the chart's measure. Also against the closed
    form b - a, pi (b^2 - a^2), 4/3 pi (b^3 - a^3) in mpmath. 16 (p + 1)^2 ulp.
    First reds: the density's factor d dropped; the mass rule scaled by half a
    panel; a measure in r instead of r^d.
    """
    geometry = _BODIES[name]
    basis = _basis(geometry, degree, layers, ratio)
    W = np.asarray(basis.mass)
    p1 = degree + 1
    ends = np.asarray(basis.partition.breakpoints)
    chart = Chart(geometry.coord)
    measure = chart.measure(ends)
    tol = 16 * p1 * p1 * _EPS
    for p in range(len(ends) - 1):
        total = W[p * p1:(p + 1) * p1, p * p1:(p + 1) * p1].sum()
        np.testing.assert_allclose(total, measure[p], rtol=tol, atol=0.0, err_msg=f"panel {p}")
        np.testing.assert_allclose(total, float(R.volume(_chart_name(name), ends[p], ends[p + 1])), rtol=tol, atol=0.0)
    np.testing.assert_allclose(W.sum(), chart.measure(np.asarray(geometry.breakpoints)).sum(), rtol=tol, atol=0.0)


# ── EB1-EB3: the even panel at a singular stratum (rung 3) ─────────────────
#
# P1 step (b), third rung (the plan's API sketch item 3 and the user's ruling of
# 2026-10-06; spec ``scratch/characteristic_architecture/p1_step_b3/spec.md``
# rows EB1-EB3): on the panel whose lower end is the centre or the axis, the
# functions are Lagrange polynomials in c^2 through the squares of the same
# Gauss points (Schwarz: a smooth O(d)-invariant function is a smooth function
# of |x|^2), and the mass rule takes 2p + 2 points on every panel.


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree", "layers", "ratio"), _ARGS, ids=_IDS)
@pytest.mark.rests_on(_HERE + "test_every_breakpoint_is_a_panel_end_and_each_panel_lies_in_one_region")
def test_the_even_panel_is_the_one_touching_a_singular_stratum(name, degree, layers, ratio) -> None:
    """[EB1] ``even`` is true exactly on panel 0 of a solid sphere or cylinder; on no panel of a hollow body or a slab.

    Written from the body (its first breakpoint and its chart), not from
    ``singular_strata``. First red: ``even`` keyed on the first breakpoint being
    0 whatever the chart (true on the slab, whose x = 0 is a face).
    """
    geometry = _BODIES[name]
    basis = _basis(geometry, degree, layers, ratio)
    want = np.zeros(basis.n_panels, dtype=bool)
    want[0] = _solid(name)
    np.testing.assert_array_equal(basis.even, want)


_EVEN_CASES = [(f"{body}-p{p}", body, p) for body in ("sphere_solid", "cylinder_solid") for p in (0, 1, 3, 5)]


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree"), [c[1:] for c in _EVEN_CASES], ids=[c[0] for c in _EVEN_CASES])
@pytest.mark.rests_on(_HERE + "test_the_even_panel_is_the_one_touching_a_singular_stratum",
                      _HERE + "test_the_functions_are_cardinal_at_the_nodes_and_sum_to_one")
def test_the_even_panel_spans_the_polynomials_in_c_squared(name, degree) -> None:
    """[EB2] On the even panel an even polynomial of degree 2p, sampled at the nodes, is reproduced; c^(2p+1) is not.

    The polynomial (seeded coefficients of 1, c^2, ..., c^2p) is read at both
    panel ends, a hair inside them and seven random points, against mpmath: 64
    ulp of its maximum (`[M]` 2026-10-06 on the spec's prototype: 0, 0.6, 5.7,
    37.6 ulp at p = 0, 1, 3, 5). The odd leg: the interpolant of c^(2p+1) (c at
    p = 0) misses it by more than 1e-3 (`[M]` 0.50, 0.19, 2.7e-2, 3.7e-3). First
    red: the panel's functions Lagrange in c (the rung-2 basis): the even
    polynomial misses (`[M]` 0.13, 9.1e-5, 1.5e-6 at p = 1, 3, 5) and the odd leg
    is reproduced.
    """
    basis = _basis(_BODIES[name], degree, 0, 0.5)
    a, b = basis.partition.breakpoints[:2]
    rng = np.random.default_rng(40 + degree)
    coefficients = rng.uniform(-1.0, 1.0, degree + 1)
    even = [0.0] * (2 * degree + 1)
    even[::2] = list(coefficients)
    nodes = np.asarray(basis.nodes[:degree + 1])
    c = a + (b - a) * np.concatenate([[0.0, 1.0, 1e-9, 1.0 - 1e-9], rng.random(7)])
    panel = np.zeros(c.shape, dtype=int)

    def value(poly, x):
        return np.array([float(mp.polyval(list(reversed([R.mpf(v) for v in poly])), R.mpf(t))) for t in x])

    got = basis.values(c, panel) @ value(even, nodes)
    want = value(even, c)
    assert np.max(np.abs(got - want)) <= 64 * _EPS * np.max(np.abs(want)), np.max(np.abs(got - want)) / np.max(np.abs(want)) / _EPS
    odd = [0.0] * (2 * degree + 1) + [1.0]
    miss = np.max(np.abs(basis.values(c, panel) @ value(odd, nodes) - value(odd, c))) / np.max(np.abs(value(odd, c)))
    assert miss > 1e-3, f"c^{2 * degree + 1} is reproduced to {miss:.1e}: the panel holds an odd mode"


_EVEN_MASS = [(f"{body}-p{p}", body, p) for body in ("sphere_solid", "cylinder_solid") for p in (1, 3, 5)]


@pytest.mark.foundation
@pytest.mark.parametrize(("name", "degree"), [c[1:] for c in _EVEN_MASS], ids=[c[0] for c in _EVEN_MASS])
@pytest.mark.rests_on(_HERE + "test_the_even_panel_spans_the_polynomials_in_c_squared",
                      _HERE + "test_the_volume_density_is_the_charts_by_hand")
def test_the_even_panels_mass_is_the_volume_integral_of_its_even_products(name, degree) -> None:
    """[EB3] The even panel's mass block against mpmath products of ``R.lagrange_even`` and the hand-written density.

    64 (1 + kappa) ulp of the block's largest entry (`[M]` 2026-10-06 on the
    spec's prototype: 5.6 to 9.3 ulp). First red: the mass rule at p + 2 points
    (the rung-2 rule): `[M]` red on the sphere at p = 1, 3, 5 and on the cylinder
    at p = 3, 5; BLIND on the cylinder at p = 1 (degree 4p + 1 = 5 is exact at
    3 points), declared, so the cylinder's catchers are p >= 3.
    """
    chart = _chart_name(name)
    basis = _basis(_BODIES[name], degree, 0, 0.5)
    p1 = degree + 1
    a, b = basis.partition.breakpoints[:2]
    got = np.asarray(basis.mass)[:p1, :p1]
    want = _mass_reference(chart, a, b, degree)
    scale = float(np.max(np.abs(want))) * (1.0 + _conditioning(a, b))
    assert np.max(np.abs(got - want)) <= 64 * _EPS * scale, np.max(np.abs(got - want)) / scale / _EPS


# ── P9: the basis reads only the partition and the resolution (C6b) ───────


@pytest.mark.foundation
def test_the_basis_reads_only_the_bodys_partition_and_the_resolution() -> None:
    """[P9, C6b] Two bodies differing only in materials and laws give one basis, bit for bit; the factory takes no points.

    The factory's parameters are exactly (regions, degree, layers, ratio): no
    reading point and no line can enter the basis (C6b, the knots do not move
    with the evaluation points). A deliberately wrong structure as the control:
    a body with one breakpoint moved gives different panel ends.
    """
    assert list(inspect.signature(PanelBasis.of).parameters) == ["regions", "degree", "layers", "ratio"]
    a = StructuredGeometry.sphere(_MR3, (0, 1, 2), outer=_spec(0.6))
    b = StructuredGeometry.sphere(_MR3, (2, 0, 7), outer=VacuumInflow())
    moved = StructuredGeometry.sphere((0.0, 0.5, 1.4, 2.0), (0, 1, 2), outer=_spec(0.6))
    ba, bb, bm = (_basis(g, 3, 2, 0.5) for g in (a, b, moved))
    assert ba.partition.breakpoints == bb.partition.breakpoints
    np.testing.assert_array_equal(ba.nodes, bb.nodes)
    np.testing.assert_array_equal(ba.mass, bb.mass)
    assert ba.partition.breakpoints != bm.partition.breakpoints


# ── the refinement's refusals ─────────────────────────────────────────────


@pytest.mark.foundation
def test_a_basis_refuses_a_partition_that_does_not_refine_the_body_with_its_ends() -> None:
    """A directly built basis refuses panel ends missing a breakpoint, with other ends, or on another chart.

    Each refusal keyed to its own fragment. First reds: the refinement check
    deleted (a panel straddling an interface is then constructible).
    """
    regions = ConcentricPartition.of(StructuredGeometry.sphere(_MR3, (0, 1, 2), outer=_spec(0.6)))
    good = PanelBasis.of(regions, 2, 1, 0.5)
    sphere = regions.chart

    def panels(ends, chart=sphere):
        return PanelBasis(regions, ConcentricPartition(chart, tuple(ends)), 2)

    assert panels(good.partition.breakpoints).partition.breakpoints == good.partition.breakpoints
    with pytest.raises(ValueError, match="refine the body's breakpoints"):
        panels((0.0, 0.3, 1.0, 1.7, 2.0))                    # 0.5 and 1.5 missing: panels straddle
    with pytest.raises(ValueError, match="refine the body's breakpoints"):
        panels((0.0, 0.5, 1.5, 2.0, 2.4))                    # another outer end
    with pytest.raises(ValueError, match="posed on one chart"):
        panels(regions.breakpoints, chart=Chart(CoordSystem.CYLINDRICAL))


# ── W1, W2: the walls re-keyed onto the panel partition ──────────────────


def _walls_cases():
    sph = StructuredGeometry.sphere(_MR3, (0, 1, 2), outer=_spec(0.6))
    hol = StructuredGeometry.sphere(_MR3H, (0, 1, 2), inner=_spec(0.3), outer=_spec(0.8))
    slab = StructuredGeometry.slab(_SLB3, (0, 1, 2), left=_spec(0.3), right=_spec(0.8))
    per = StructuredGeometry.slab(_SLB3, (0, 1, 2), left=PeriodicBoundary(axis="x"), right=PeriodicBoundary(axis="x"))
    # (id, geometry, [(breakpoint on the panels' index: "first" / "last", specular, partner: "self" / "other")])
    return [
        ("sphere_solid", sph, [("last", 0.6, "self")]),
        ("sphere_hollow", hol, [("first", 0.3, "self"), ("last", 0.8, "self")]),
        ("slab_mirrors", slab, [("first", 0.3, "self"), ("last", 0.8, "self")]),
        ("slab_periodic", per, [("first", 1.0, "other"), ("last", 1.0, "other")]),
    ]


_WALLS_CASES = _walls_cases()


@pytest.mark.foundation
@pytest.mark.parametrize(("geometry", "expected"), [c[1:] for c in _WALLS_CASES], ids=[c[0] for c in _WALLS_CASES])
@pytest.mark.rests_on(_HERE + "test_every_breakpoint_is_a_panel_end_and_each_panel_lies_in_one_region",
                      _WALLS + "test_every_law_reads_as_the_wall_its_physics_names")
def test_the_walls_rekey_onto_the_panel_partition(geometry, expected) -> None:
    """[W1] Walls.on(panels): breakpoint 0 stays 0, n becomes the panel count P, partners follow, amplitudes kept.

    The expected walls are written from the body (which end, which amplitude,
    whether a wrap pairs the two ends), never from ``Walls.of``'s output. The
    re-keyed walls answer the per-breakpoint lookups at 0 and P and refuse
    every interior panel end (IndexError, as an interior breakpoint does). First
    reds: n kept as the last index (the outer wall lost); a wrap's partner left
    at n (the periodic slab's cycle broken); the amplitudes of the two walls
    swapped.
    """
    basis = _basis(geometry, 3, 2, 0.5)
    P = len(basis.partition.breakpoints) - 1
    walls = Walls.of(geometry).on(basis.partition)
    assert walls.n_regions == P
    index = {"first": 0, "last": P}
    got = sorted((w.breakpoint, w.specular, w.diffuse, w.partner) for w in walls.walls)
    want = sorted((index[end], amp, 0.0, index[end] if partner == "self" else P - index[end])
                  for end, amp, partner in expected)
    assert got == want
    wall_points = np.array([b for b, *_ in want])
    np.testing.assert_array_equal(walls.specular_at(wall_points), [a for _, a, *_ in want])
    for interior in (1, P // 2, P - 1):
        with pytest.raises(IndexError):
            walls.specular_at(np.array(interior))


@pytest.mark.foundation
def test_walls_refuse_a_partition_on_another_chart() -> None:
    """[W2] Walls.on refuses a partition posed on another chart, naming both charts.

    Declared gap (spec, finding F1): ``Walls`` holds no breakpoints, so it
    cannot refuse a partition with other ends or missing a breakpoint; that
    refusal is the panel basis's (``test_a_basis_refuses_a_partition_that_does_not_refine_the_body_with_its_ends``).
    """
    geometry = StructuredGeometry.sphere(_MR3, (0, 1, 2), outer=_spec(0.6))
    cylinder = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), _MR3)
    with pytest.raises(ValueError, match="are not re-keyed onto a"):
        Walls.of(geometry).on(cylinder)
