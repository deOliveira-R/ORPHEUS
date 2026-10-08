"""Gates for the characteristic reference's door (:class:`~orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation`).

P1 step (b), fifth rung, first half (5a), of the characteristic-reference
campaign (``.claude/plans/characteristic_reference_architecture.md``, "P1 step
(b), fifth rung: API sketch", ruled 2026-10-07: the rung splits; a ``Response``
is the group-transposed forward problem with the detector's retraction as its
source, so its answer is the adjoint scalar flux; the eigen flux of a finite body
is gauged by its total production, sum_g int nu Sigma_f phi dV = 1 (the rulings of
the same day, after qa's review, withdrawing the density gauge 100);
``Nearest(tau)`` is served with tau in the k chart and reads its eigenvalue only). Verification spec
``scratch/characteristic_architecture/p1_step_b5/spec.md``: rows P1-P4, D1-D10,
K1-K4, N1-N4, F1-F5, A1-A2, R1-R3, C1-C5.

The ladder, bottom up (``rests_on`` on each row):

1. rung 4, landed: the Galerkin system's pencil, source questions and adjoint
   (``test_characteristic_system.py``);
2. the projection of a function onto the panel basis [P1-P4];
3. the door's refusals and its laziness [D1-D10];
4. the reading of a flux integral against the Galerkin flux [R1, R2];
5. the questions through the door: the eigenvalue [K1-K4], the mode nearest
   tau [N1-N4], the fixed source [F1-F5], the response [A1, A2];
6. the reference solution and the traced memo [R3, C1-C5].

Every closed form is written here from the ``Mixture`` arrays or the geometry's
numbers (volumes, the mirror slab's modes), never through the door's objects
(``instrument-doctrine`` X4); the rows that read the door's own flux
coefficients (R1) say so, since their claim is the reading and not the solve.
Every value row reads in this process (:func:`~orpheus.numerics.traced_memo.bypass`),
so one solve serves a derivation's readings; the memo's own rows (C2-C5) run
through it. Bands are measured (`[M]` 2026-10-07, ``.venv/bin/python -O``,
``scratch/characteristic_architecture/p1_step_b5/ta/``); each row's first red is a
battery arm (``scratch/characteristic_architecture/p1_step_b5/ta/battery/``).
"""
from __future__ import annotations

from collections.abc import Callable
from functools import lru_cache
from typing import Any

import numpy as np
import pytest
import sympy

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.common.dense_pencil import NoLeastSolution
from orpheus.derivations.common.exact_homogeneous import exact_infinite_medium_reference
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.derivations.continuous.characteristic import (
    CharacteristicDerivation, GalerkinSystem, PanelBasis, RegionCrossSections, Resolution, TransportResolution,
    characteristic_reference,
)
from orpheus.derivations.continuous.characteristic import reference as reference_module
from orpheus.derivations.continuous.characteristic.assembly import LineRule
from orpheus.derivations.continuous.sood_registry import sood2003 as sood
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.numerics.mesh_free_function import MeshFreeFunction, RegionwiseConstant, Symbolic
from orpheus.numerics.observable import Eigenvalue, FluxIntegral, PointValue, Ratio
from orpheus.numerics.question import Eigen, FixedSource, Nearest, Response
from orpheus.numerics.traced_memo import bypass
from orpheus.reference.reading import Uncertified
from orpheus.reference.solution import ReferenceSolution
from orpheus.specification.specification import GeometrySpecification, InfiniteMediumSpecification
from tests.gates._content_identity_helpers import (
    Entry, check_equal_pair, check_perturbation, check_population, leg, pair_ids, param_id, perturbation_ids,
)
from tests.gates.derivations.test_characteristic_system import (
    _ABS, _ADJ_ABSORBER, _ADJ_FUEL, _ADJ_REFLECTOR, _MIRROR, _N2N, _PU2, _SELF_ONLY, _UP2N, _URRB, _VACUUM,
    _fission_by_hand, _scattering_by_hand, _walls, _zero_d,
)
from tests.gates.numerics import _traced_memo_api as memo_api

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_HERE = "tests/gates/derivations/test_characteristic_reference.py::"
_SYSTEM = "tests/gates/derivations/test_characteristic_system.py::"
_D5 = _SYSTEM + "test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio"
_E2 = _SYSTEM + "test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux"
_E6 = _SYSTEM + "test_a_source_in_a_body_supercritical_by_its_n2n_emission_is_refused"
_D13IV = _SYSTEM + "test_the_adjoint_is_the_forward_solve_of_the_group_transposed_problem"
_D2 = _SYSTEM + "test_k_is_one_at_soods_two_group_critical_sizes"
_SU2 = _SYSTEM + "test_a_field_off_the_support_is_refused_and_answered_once_the_support_is_widened"
_D7 = _SYSTEM + "test_the_fundamental_mode_satisfies_its_pencil"

_K = CellCoefficient.every(Channel.FISSION_EMISSION)
#: The working point of rung 4 (``_WORK``: degree 3, 2 layers of ratio 0.4, 8 line points, 12 points along each line),
#: with 8 = 2(p + 1) projection points, the mass matrix's own rule.
_RES = Resolution(3, 2, 0.4, TransportResolution(8, 12, 12), 8)
#: The smallest resolution a body solves at, for the rows that count interpreters and read no value.
_TINY = Resolution(1, 0, 0.4, TransportResolution(4, 4, 4), 4)
_COORD = {"sphere": CoordSystem.SPHERICAL, "cylinder": CoordSystem.CYLINDRICAL, "slab": CoordSystem.CARTESIAN}
_MR3 = (0.0, 0.5, 1.5, 2.0)
_MR3H = (0.4, 0.5, 1.5, 2.0)
_SLB3 = (0.0, 0.4, 1.5, 2.3)
_ADJ_BREAKPOINTS = (0.0, 0.7, 1.1, 1.8)


# ── the bodies, as specifications ────────────────────────────────────────


def _sphere(breakpoints, mat_ids, outer: BC, inner: BC | None = None) -> StructuredGeometry:
    return StructuredGeometry.sphere(breakpoints, mat_ids, outer=outer, inner=inner)


def _slab(breakpoints, mat_ids, left: BC, right: BC) -> StructuredGeometry:
    return StructuredGeometry.slab(breakpoints, mat_ids, left=left, right=right)


def _spec(geometry: StructuredGeometry, mixtures: dict[int, Any], question: Any) -> GeometrySpecification:
    return GeometrySpecification(Materials(mixtures), geometry, question)


def _layered(geometry_of: Callable[..., StructuredGeometry], breakpoints, *laws: BC) -> StructuredGeometry:
    """One mixture under one material id per interval (a layered body, the multi-region route)."""
    return geometry_of(breakpoints, tuple(range(len(breakpoints) - 1)), *laws)


def _same(mixture, n: int) -> dict[int, Any]:
    return {i: mixture for i in range(n)}


@lru_cache(maxsize=None)
def _derivation(specification: GeometrySpecification, resolution: Resolution = _RES) -> CharacteristicDerivation:
    return CharacteristicDerivation(specification, resolution)


def _read(specification: GeometrySpecification, observable: Any, resolution: Resolution = _RES) -> float:
    """The door's reading of ``observable``, in this process."""
    with bypass():
        reading = _derivation(specification, resolution).evaluate(observable)
    assert type(reading) is Uncertified
    return float(reading.value)


def _volumes(chart: str, breakpoints) -> np.ndarray:
    """Each interval's measure, by hand: a sphere's shell, a cylinder's annulus per unit height, a slab's width per unit area."""
    r = np.asarray(breakpoints, dtype=float)
    if chart == "sphere":
        return 4.0 * np.pi * (r[1:] ** 3 - r[:-1] ** 3) / 3.0
    if chart == "cylinder":
        return np.pi * (r[1:] ** 2 - r[:-1] ** 2)
    return r[1:] - r[:-1]


def _indicator(n_regions: int, region: int, group: int, groups: int = 2) -> FluxIntegral:
    values = np.zeros((n_regions, groups))
    values[region, group] = 1.0
    return FluxIntegral(RegionwiseConstant(values))


def _rwc(values) -> RegionwiseConstant:
    return RegionwiseConstant(np.asarray(values, dtype=float))


def _symbolic_constant(values) -> Symbolic:
    return Symbolic.of(*(sympy.Float(v) if v else sympy.Integer(0) for v in values))


# ── 5a.1 the projection onto the panel basis ─────────────────────────────


@lru_cache(maxsize=None)
def _basis(chart: str, breakpoints, degree: int = 3, layers: int = 2, ratio: float = 0.4) -> PanelBasis:
    return PanelBasis.of(ConcentricPartition(Chart(_COORD[chart]), tuple(breakpoints)), degree, layers, ratio)


def _panel_of(basis: PanelBasis, c: np.ndarray) -> np.ndarray:
    ends = np.asarray(basis.partition.breakpoints)
    return np.clip(np.searchsorted(ends, c, side="right") - 1, 0, basis.n_panels - 1)


def _reconstruct(basis: PanelBasis, coefficients: np.ndarray, c: np.ndarray) -> np.ndarray:
    """sum_i coefficients[..., i] u_i(c), ``(..., q)``: the basis function's values written through ``values``."""
    panel = _panel_of(basis, c)
    u = basis.values(c, panel)                                           # (q, p + 1)
    picked = coefficients[..., basis.columns(panel)]                      # (..., q, p + 1)
    return np.einsum("...qi,qi->...q", picked, u)


_PROJECTION_BODIES = [("sphere", _MR3), ("slab", _SLB3), ("sphere", _MR3H), ("cylinder", (0.0, 0.6, 1.0))]
_PROJECTION_IDS = ["sphere3", "slab3", "hollow-sphere3", "cylinder2"]


@pytest.mark.verifies("characteristic-projection")
@pytest.mark.l0
@pytest.mark.parametrize(("chart", "breakpoints"), _PROJECTION_BODIES, ids=_PROJECTION_IDS)
def test_a_per_region_constant_projects_to_its_value_at_every_node(chart: str, breakpoints) -> None:
    """[P1; l0] ``project`` of a per-region constant (two groups, distinct values per region) returns the
    constant at every node of its region, to 1e-14 relative.

    The constant is in every panel's space (the even panel's included), so
    the L2 projection is the identity on it when the load is integrated as
    exactly as the mass matrix (2(p + 1) points). `[M]` 2026-10-07: 3.3e-15 on the
    solid sphere, 6.7e-16 on the other three (``ta/measure_1.log``). First reds: the load returned without W^-1; the
    measure density dropped from the load (the sphere's and cylinder's
    projections move O(1), the slab's measure is 1: the slab row is blind to it).
    """
    basis = _basis(chart, breakpoints)
    table = np.array([[1.0, 0.3], [2.5, 0.7], [0.2, 4.0]])[: len(breakpoints) - 1]
    interior = np.asarray(breakpoints[1:-1])
    projected = basis.project(lambda c: table[np.searchsorted(interior, c, side="right")].T, 2 * basis.per_panel)
    expected = table[basis.region].T
    assert projected.shape == (2, basis.size)
    assert np.max(np.abs(projected / expected - 1.0)) < 1e-14


@pytest.mark.verifies("characteristic-projection")
@pytest.mark.l0
@pytest.mark.parametrize(("chart", "breakpoints"), _PROJECTION_BODIES, ids=_PROJECTION_IDS)
@pytest.mark.rests_on(_HERE + "test_a_per_region_constant_projects_to_its_value_at_every_node")
def test_a_function_of_the_panel_space_projects_to_its_coefficients(chart: str, breakpoints) -> None:
    """[P2; l0] A random element of the panel space (seeded coefficients, two groups), written through
    ``values``, projects back to its coefficients to 2e-14 relative to the largest; and ``f = c``, a polynomial
    of degree 1, is reproduced at every node except on the even panel, whose space is the polynomials in c^2.

    The even leg is a declared NON-exactness, bounded one-sidedly: on the
    sphere's and the cylinder's centre panel the nodal error of ``f = c`` is
    above 1e-4 (`[M]` 2026-10-07: 2.8e-2, 2.4e-2); every other panel of every
    body reproduces it to 2e-14 of the largest node (`[M]` <= 4.5e-16; the
    coefficients <= 1.9e-15). First reds: the load integrated with p + 1
    points (the even panel's products are degree 4p + d - 1); the even panel
    read as an ordinary polynomial panel (``even`` all False: the c leg then
    reproduces c on the centre panel and its bound reddens).
    """
    basis = _basis(chart, breakpoints)
    points = 2 * basis.per_panel
    coefficients = np.random.default_rng(31).random((2, basis.size))
    recovered = basis.project(lambda c: _reconstruct(basis, coefficients, c), points)
    assert np.max(np.abs(recovered - coefficients)) < 2e-14 * np.max(np.abs(coefficients))
    linear = basis.project(lambda c: c[None, :], points)[0]
    error = np.abs(linear - basis.nodes)
    even_nodes = basis.even[basis.panel]
    assert np.max(error[~even_nodes]) < 2e-14 * np.max(basis.nodes)
    if even_nodes.any():
        assert np.max(error[even_nodes]) > 1e-4, "the centre panel reproduced an odd function: is it still even?"
    else:
        assert chart == "slab" or breakpoints[0] > 0.0                   # only a solid round body has an even panel


def _smooth(c: np.ndarray) -> np.ndarray:
    """A smooth non-polynomial function of two groups, even in c (so smooth at a solid body's centre)."""
    return np.stack([np.exp(-c**2), 1.0 / (1.0 + c**2)])


@pytest.mark.verifies("characteristic-projection")
@pytest.mark.l0
@pytest.mark.parametrize(("chart", "breakpoints"), _PROJECTION_BODIES[:2], ids=_PROJECTION_IDS[:2])
@pytest.mark.rests_on(_HERE + "test_a_function_of_the_panel_space_projects_to_its_coefficients")
def test_the_projection_of_a_smooth_function_converges_with_the_degree(chart: str, breakpoints) -> None:
    """[P3; l0] The reconstruction error of a smooth non-polynomial function (max over 400 points per panel)
    falls by a factor above 3 with every degree from 1 to 6, and is below 1e-6 at degree 6 (the panels of the
    working grading: 2 layers of ratio 0.4).

    `[M]` 2026-10-07: sphere 2.3e-2 to 2.9e-7, slab 2.5e-2 to 4.6e-7, each step
    a factor 5.0 or more. First red: the load integrated with a fixed 2-point
    rule (the error stalls).
    """
    errors = []
    for degree in range(1, 7):
        basis = _basis(chart, breakpoints, degree)
        coefficients = basis.project(_smooth, 2 * basis.per_panel)
        ends = np.asarray(basis.partition.breakpoints)
        c = np.concatenate([np.linspace(a, b, 402)[1:-1] for a, b in zip(ends[:-1], ends[1:])])
        errors.append(float(np.max(np.abs(_reconstruct(basis, coefficients, c) - _smooth(c)))))
    assert all(b < a / 3.0 for a, b in zip(errors[:-1], errors[1:])), errors
    assert errors[-1] < 1e-6, errors


@pytest.mark.verifies("characteristic-projection")
@pytest.mark.l0
@pytest.mark.rests_on(_HERE + "test_a_function_of_the_panel_space_projects_to_its_coefficients")
def test_a_step_inside_a_panel_is_integrated_once_it_is_named() -> None:
    """[P4; l0] The indicator of c < s, s inside a panel of the sphere's middle region: ``project`` with
    ``steps=(s,)`` equals W_P^-1 of the load integrated here piece by piece (numpy's Gauss-Legendre, 20 points,
    the panel split at s) to 1e-13; without the step it differs by more than 1e-3 (the step matters there).

    `[M]` 2026-10-07: 4.0e-15 with the step, 0.22 without. First red: the
    steps ignored.
    """
    basis = _basis("sphere", _MR3)
    ends = np.asarray(basis.partition.breakpoints)
    split = int(np.searchsorted(ends, 1.0))                               # a panel in the middle region
    s = 0.5 * (ends[split - 1] + ends[split]) + 0.1 * (ends[split] - ends[split - 1])
    target = split - 1
    nodes, weights = np.polynomial.legendre.leggauss(20)

    def load_on(a: float, b: float, panel: int) -> np.ndarray:
        c = 0.5 * (b - a) * nodes + 0.5 * (a + b)
        return basis.values(c, np.full(c.shape, panel)).T @ (0.5 * (b - a) * weights * 4.0 * np.pi * c**2)

    expected = np.zeros(basis.size)
    for panel in range(basis.n_panels):
        a, b = ends[panel], ends[panel + 1]
        if b <= s:
            load = load_on(a, b, panel)
        elif panel == target:
            load = load_on(a, s, panel)
        else:
            continue
        columns = basis.columns(np.array(panel))
        expected[columns] = np.linalg.solve(basis.mass[np.ix_(columns, columns)], load)

    def step(c: np.ndarray) -> np.ndarray:
        return (c < s).astype(float)[None, :]

    named = basis.project(step, 8, (s,))[0]
    unnamed = basis.project(step, 8)[0]
    assert np.max(np.abs(named - expected)) < 1e-13
    assert np.max(np.abs(unnamed - expected)) > 1e-3


@pytest.mark.l0
def test_the_cross_sections_carry_chi_and_production_and_transpose_by_swapping_them() -> None:
    """[X1; l0] ``production`` is nu Sigma_f (``SigP``) region by region, bit for bit; ``transposed()`` keeps the
    total, transposes the scattering bit for bit, and swaps the fission's two factors: its ``production`` is
    chi and its fission is the forward fission transposed. On ``_UP2N`` | ``_ABS`` | ``_PU2`` (chi in both
    groups, a non-fissile region).

    Fix (2) of 2026-10-07 stores chi and nu Sigma_f and derives the
    fission. `[M]` before the fix: ``production`` was ``fission.sum(axis=1)``
    (nu Sigma_f times sum chi, 1 ulp off) and the transposed production chi
    times sum nu Sigma_f; both legs were red. First reds: the production read
    as ``fission.sum(axis=2)`` (chi times sum nu Sigma_f); ``transposed()``
    leaving the fission untransposed.
    """
    mixtures = [_UP2N, _ABS, _PU2]
    xs = RegionCrossSections.of(mixtures)
    adjoint = xs.transposed()
    for r, m in enumerate(mixtures):
        np.testing.assert_array_equal(xs.production[r], m.SigP)
        np.testing.assert_array_equal(adjoint.production[r], m.chi)
        np.testing.assert_array_equal(adjoint.total[r], xs.total[r])
        np.testing.assert_array_equal(adjoint.scattering[r], xs.scattering[r].T)
        np.testing.assert_array_equal(adjoint.fission[r], xs.fission[r].T)


# ── 5a.2 the door: refusals and laziness ─────────────────────────────────


def _hetero_sphere(question: Any) -> GeometrySpecification:
    """A three-region vacuum sphere, material ids out of interval order, both materials fissile."""
    return _spec(_sphere(_MR3, (1, 0, 1), BC.vacuum), {0: _ABS, 1: _UP2N}, question)


def _two_fissile(question: Any) -> GeometrySpecification:
    return _spec(_sphere(_MR3, (1, 0, 1), BC.vacuum), {0: _PU2, 1: _UP2N}, question)


def _count_line_rules(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """A spy on ``LineRule.of``, the start of every transport block's assembly: one entry per call."""
    calls: list[int] = []
    original: Any = LineRule.of

    def counting(cls: Any, *args: Any, **kwargs: Any) -> Any:
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(LineRule, "of", classmethod(counting))
    return calls


_MU_SOURCE = Symbolic.of(1 + Symbolic.mu, 0)
_POINT = {CellCoefficient.every(Channel.SCATTERING_EMISSION): 0.0}

_REFUSALS: list[tuple[str, Callable[[], Any], type, str]] = [
    ("infinite-medium", lambda: InfiniteMediumSpecification(0, _UP2N, Eigen(_K)), TypeError,
     "an infinite medium has no geometry"),
    ("eigen-along-scattering", lambda: _hetero_sphere(Eigen(CellCoefficient.every(Channel.SCATTERING_EMISSION))),
     ValueError, "answers an eigenvalue along the fission emission"),
    ("eigen-along-one-fissile-material", lambda: _two_fissile(Eigen(CellCoefficient(((0, Channel.FISSION_EMISSION),)))),
     ValueError, "answers an eigenvalue along the fission emission"),
    ("eigen-off-the-physical-point", lambda: _hetero_sphere(Eigen(_K, point=_POINT)), ValueError,
     "answers at the physical point"),
    ("source-off-the-physical-point", lambda: _hetero_sphere(FixedSource(_rwc(np.ones((3, 2))), point=_POINT)),
     ValueError, "answers at the physical point"),
    ("response-off-the-physical-point", lambda: _hetero_sphere(Response(_rwc(np.ones((3, 2))), point=_POINT)),
     ValueError, "answers at the physical point"),
    ("anisotropic-scattering-at-interval-2", lambda: _spec(
        _sphere(_MR3, (1, 1, 0), BC.vacuum), {0: get_mixture("B", "2g"), 1: _PU2}, Eigen(_K)),
     NotImplementedError, "region 2's SigS has a non-zero Legendre order 1"),
    ("source-depending-on-mu-sphere", lambda: _hetero_sphere(FixedSource(_MU_SOURCE)), NotImplementedError,
     "serves an isotropic source only"),
    ("source-depending-on-phi-slab", lambda: _spec(
        _slab(_SLB3, (1, 0, 1), BC.reflective, BC.vacuum), {0: _ABS, 1: _UP2N},
        FixedSource(Symbolic.of(1 + Symbolic.phi**2, 0))), NotImplementedError, "serves an isotropic source only"),
    ("detector-depending-on-mu", lambda: _hetero_sphere(Response(_MU_SOURCE)), NotImplementedError,
     "serves an isotropic detector only"),
    ("unknown-wall-tag", lambda: _spec(_sphere(_MR3, (1, 0, 1), BC("marshak")), {0: _ABS, 1: _UP2N}, Eigen(_K)),
     NotImplementedError, "the reference admits the tag kinds"),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("build", "error", "fragment"), [r[1:] for r in _REFUSALS], ids=[r[0] for r in _REFUSALS])
def test_the_door_refuses_at_construction_naming_what_it_refused(build, error: type, fragment: str,
                                                                  monkeypatch: pytest.MonkeyPatch) -> None:
    """[D1-D7] Each refusal is raised by the constructor, with its exception type and its shortest distinctive
    fragment, before any transport block is assembled (a spy on ``LineRule.of`` counts 0).

    The fragments are pairwise disjoint over the rows that differ in
    mechanism (asserted below), so a refusal reached by another guard reddens
    its row. The anisotropy row puts the anisotropic material at INTERVAL 2
    under material id 0: the fragment names the interval, so a door reading
    the mixtures in id order (region 0) reddens it. First reds, one per arm:
    each guard removed; the mixtures read by id order; the point guard keyed
    on ``Eigen`` only; the direction guard on the source only; the walls read
    lazily (the tag row then constructs).
    """
    calls = _count_line_rules(monkeypatch)
    specification = build()
    with pytest.raises(error, match=fragment):
        CharacteristicDerivation(specification, _RES)
    assert calls == []


def test_the_refusal_fragments_are_disjoint() -> None:
    """[D1-D7, the instrument] No two mechanisms share a fragment: a fragment matching another row's
    mechanism would let a wrong guard satisfy the row."""
    mechanisms = {fragment for *_, fragment in _REFUSALS}
    for a in mechanisms:
        for b in mechanisms - {a}:
            assert a not in b


@pytest.mark.foundation
@pytest.mark.parametrize("question", [Eigen(_K), Eigen(_K, mode=Nearest(0.5)), FixedSource(_rwc(np.ones((3, 2)))),
                                      Response(_rwc(np.ones((3, 2))))], ids=["eigen", "nearest", "source", "response"])
def test_a_served_question_constructs_and_assembles_nothing(question, monkeypatch: pytest.MonkeyPatch) -> None:
    """[D8] The positive leg of D1-D7: every served question constructs (the k question spelled UNRESOLVED,
    as the specification canonicalises it), and construction assembles no transport block (``LineRule.of``
    0 calls); the first reading assembles them (the activation leg: the count is then positive).

    First red: the system's blocks built in the constructor.
    """
    calls = _count_line_rules(monkeypatch)
    spec = _hetero_sphere(question)
    derivation = CharacteristicDerivation(spec, _TINY)
    assert calls == []
    with bypass():
        derivation.evaluate(Eigenvalue() if isinstance(question, Eigen) else _indicator(3, 0, 0))
    assert calls, "the activation leg: the reading assembled nothing"


@pytest.mark.foundation
@pytest.mark.parametrize(("fields", "error", "fragment"), [
    ({"transport": (8, 12, 12)}, TypeError, "is a TransportResolution"),
    ({"source_points": 7}, ValueError, "takes at least"),
    ({"source_points": 0}, ValueError, "takes at least"),
    ({"degree": 3.0}, TypeError, "is an integer"),
    ({"layers": 2.0}, TypeError, "is an integer"),
    ({"source_points": 8.0}, TypeError, "is an integer"),
], ids=["untyped-transport", "fewer-points-than-the-mass-rule", "no-projection-point", "real-degree", "real-layers",
        "real-source-points"])
def test_a_malformed_resolution_is_refused(fields: dict, error: type, fragment: str) -> None:
    """[D9] ``Resolution`` refuses a transport resolution that is not one, fewer projection points than the
    mass matrix's rule, 2(p + 1) = 8 at degree 3 (below it a symbolic function of the panel space is no longer
    projected to rounding, P2), and a degree, a number of layers or of points that is not an integer.
    First reds: each guard removed."""
    values = {"degree": 3, "layers": 2, "ratio": 0.4, "transport": TransportResolution(8, 12, 12), "source_points": 8}
    with pytest.raises(error, match=fragment):
        Resolution(**(values | fields))


@pytest.mark.foundation
def test_a_point_value_is_read_from_the_solved_emission_and_refused_on_a_mode_after_the_solve(monkeypatch: pytest.MonkeyPatch) -> None:
    """[D10, re-posed at rung 5b] Until 5b this row pinned ``PointValue``'s refusal naming rung 5b before any solve;
    the refusal retired with ``_refuse_point_reading`` (``test_characteristic_reading.py`` carries the reading's
    values). What survives of the door's claim: a ``PointValue`` of the fundamental is the system's reading of the
    answer's emission, bit for bit, solved once (a spy on the k pencil's fundamental counts the solve); a
    ``PointValue`` of a ``Nearest`` answer is refused naming the missing flux scale, AFTER its solve (the mode must
    be found before it is known to be a higher one).

    First reds: the reading taken from the Galerkin flux coefficients in place of the emission (the value leg);
    the refusal placed before the solve (the spy leg).
    """
    calls = []
    original = reference_module.GalerkinSystem.pencil

    def counting(self):
        calls.append(1)
        return original.func(self)

    monkeypatch.setattr(reference_module.GalerkinSystem, "pencil", property(counting))
    derivation = CharacteristicDerivation(_hetero_sphere(Eigen(_K)), _TINY)
    with bypass():
        value = derivation.evaluate(PointValue(1.0, 0)).value
    answer = derivation.answer
    assert isinstance(answer, reference_module._FundamentalAnswer)
    assert value == float(derivation.system.point_flux(1.0, answer.emission)[0])
    assert calls, "the activation leg: the point value solved nothing"
    calls.clear()
    nearest = CharacteristicDerivation(_hetero_sphere(Eigen(_K, mode=Nearest(0.5))), _TINY)
    with bypass(), pytest.raises(NotImplementedError, match="has no flux scale"):
        nearest.evaluate(PointValue(1.0, 0))
    assert calls, "the refusal came before the mode was solved"


@pytest.mark.foundation
def test_an_eigenvalue_of_a_source_answer_is_refused_before_any_solve(monkeypatch: pytest.MonkeyPatch) -> None:
    """[D11] ``evaluate(Eigenvalue())`` on a ``FixedSource`` derivation raises ValueError naming the eigen
    question, and assembles no transport block (fix (3) of 2026-10-07: the answer typed per
    question). The admission (``admit_observable``) refuses it first through ``ReferenceSolution.read``; this
    row is the derivation's own door. `[M]` red before fix (3) (the answer was solved, then the eigenvalue
    refused). First red: the refusal placed after the solve."""
    calls = _count_line_rules(monkeypatch)
    derivation = CharacteristicDerivation(_hetero_sphere(FixedSource(_rwc(np.ones((3, 2))))), _TINY)
    with bypass(), pytest.raises(ValueError, match="an eigenvalue is read from an eigen question"):
        derivation.evaluate(Eigenvalue())
    assert calls == []


# ── 5a.3 the reading of a flux integral ──────────────────────────────────


def _volume_integral(basis: PanelBasis, flux: np.ndarray, weight: Callable[[np.ndarray], np.ndarray],
                     cut: float | None = None) -> float:
    """sum_g int w_g phi_h,g dV by numpy Gauss-Legendre (24 points per panel, split at ``cut``), the measure written by hand."""
    ends = list(np.asarray(basis.partition.breakpoints))
    if cut is not None:
        ends = sorted(set(ends) | {cut})
    nodes, weights = np.polynomial.legendre.leggauss(24)
    total = 0.0
    density = {CoordSystem.SPHERICAL: lambda c: 4.0 * np.pi * c**2, CoordSystem.CYLINDRICAL: lambda c: 2.0 * np.pi * c,
               CoordSystem.CARTESIAN: lambda c: np.ones_like(c)}[basis.regions.chart.coord]
    for a, b in zip(ends[:-1], ends[1:]):
        c = 0.5 * (b - a) * nodes + 0.5 * (a + b)
        total += float(np.sum(_reconstruct(basis, flux, c) * weight(c) * (0.5 * (b - a) * weights * density(c))))
    return total


_STEP_AT = 1.13                                                            # inside a panel of the middle region of _MR3


def _weights() -> list[tuple[str, MeshFreeFunction, Callable[[np.ndarray], np.ndarray], float | None]]:
    table = np.array([[0.3, 1.1], [2.0, 0.4], [0.7, 0.9]])
    interior = np.asarray(_MR3[1:-1])
    r = Symbolic.r
    smooth = Symbolic.of(1 + r**2 * sympy.exp(-r), sympy.cos(r))
    step = sympy.Piecewise((1, r < sympy.Rational(str(_STEP_AT))), (0, True))
    return [
        ("regionwise", _rwc(table), lambda c: table[np.searchsorted(interior, c, side="right")].T, None),
        ("smooth-symbolic", smooth, lambda c: np.stack([1 + c**2 * np.exp(-c), np.cos(c)]), None),
        ("step-symbolic", Symbolic.of(step, 2 * step), lambda c: np.stack([c < _STEP_AT, 2.0 * (c < _STEP_AT)]).astype(float),
         _STEP_AT),
    ]


_SUBCRITICAL_SPHERE = _spec(_sphere(_MR3, (1, 0, 1), BC.vacuum), {0: _ABS, 1: _UP2N}, FixedSource(_rwc([[1.0, 0.4], [0.0, 0.0], [1.0, 0.4]])))


@pytest.mark.l1
@pytest.mark.parametrize("case", range(3), ids=["regionwise", "smooth-symbolic", "step-symbolic"])
@pytest.mark.rests_on(_HERE + "test_a_per_region_constant_projects_to_its_value_at_every_node",
                      _HERE + "test_a_step_inside_a_panel_is_integrated_once_it_is_named")
@pytest.mark.verifies("characteristic-door-pairing")
def test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux(case: int) -> None:
    """[R1; l1] ``FluxIntegral(w)`` equals sum_g int w_g phi_h,g dV, phi_h the door's Galerkin flux coefficients
    reconstructed through the basis and integrated here (numpy Gauss-Legendre, 24 points per panel, split at the
    weight's step), to 1e-13 relative, for a region-wise weight, a smooth symbolic weight and a symbolic step
    inside a panel.

    The pairing (W c_w)^T phi_h is the volume integral of w against phi_h
    whatever w is, when the load is integrated exactly: the reading carries
    the Galerkin flux's error, never the weight's projection error (a
    correction to the sketch's "carries the weight's projection error"). This
    row reads the door's own flux coefficients: its claim is the reading, not
    the solve. `[M]` 2026-10-07: 2.2e-16, 1.1e-16, 1.1e-16 (the smooth weight's load at 8
    points per piece is below rounding on these panels). First reds: the pairing without W (c_w^T phi); the weight
    read as a rate (4 pi); the symbolic weight's steps dropped (the step row
    moves 1e-3); the weight's groups paired with the reversed flux groups.
    """
    name, weight, by_hand, cut = _weights()[case]
    derivation = _derivation(_SUBCRITICAL_SPHERE)
    reading = _read(_SUBCRITICAL_SPHERE, FluxIntegral(weight))
    answer = derivation.answer
    assert isinstance(answer, reference_module._SourceAnswer)
    expected = _volume_integral(derivation.basis, answer.flux, by_hand, cut)
    assert abs(reading / expected - 1.0) < 1e-13, (name, reading, expected)


# ── 5a.4 the eigenvalue ──────────────────────────────────────────────────


_WIRING = [
    ("slab-mirror-vacuum", "slab", _SLB3, (_MIRROR, _VACUUM), lambda: _slab(_SLB3, (1, 0, 2), BC.reflective, BC.vacuum)),
    ("sphere-vacuum", "sphere", _MR3, (_VACUUM,), lambda: _sphere(_MR3, (1, 0, 2), BC.vacuum)),
]
_WIRING_MATERIALS = {0: _PU2, 1: _UP2N, 2: _ABS}


@pytest.mark.l1
@pytest.mark.parametrize(("chart", "breakpoints", "laws", "geometry"), [w[1:] for w in _WIRING], ids=[w[0] for w in _WIRING])
@pytest.mark.rests_on(_D7)
def test_the_doors_k_is_the_fundamental_of_the_system_written_by_hand(chart, breakpoints, laws, geometry) -> None:
    """[K1; l1] The door's ``Eigenvalue()`` is the fundamental k of the ``GalerkinSystem`` built here from the
    same numbers (the walls written as ``Wall`` tuples, the mixtures in interval order, ``PanelBasis.of`` at
    the resolution's degree and grading), bit for bit; and ``Eigenvalue()`` of a ``Nearest`` posed at that k
    reads it to 1e-13 relative.

    The material ids run (1, 0, 2) over the intervals, and the slab has its
    mirror on the LEFT. A wiring row: both sides run one assembly, so it pins
    the door's reading of the specification, not the physics (K2-K4 do).
    `[M]` 2026-10-07: equal bitwise; the Nearest leg 0. First reds:
    the mixtures read in id order; the slab's walls swapped; ``degree`` and
    ``layers`` swapped in the reading of the resolution; the eigenvalue read
    in the 1/k chart.
    """
    spec = _spec(geometry(), _WIRING_MATERIALS, Eigen(_K))
    system = GalerkinSystem(
        _basis(chart, breakpoints), _walls(chart, breakpoints, laws),
        RegionCrossSections.of([_UP2N, _PU2, _ABS]), TransportResolution(8, 12, 12),
    )
    k = system.pencil.fundamental().k
    assert _read(spec, Eigenvalue()) == k
    nearest = _spec(geometry(), _WIRING_MATERIALS, Eigen(_K, mode=Nearest(k)))
    assert abs(_read(nearest, Eigenvalue()) / k - 1.0) < 1e-13


_CLOSED_BODIES = [
    ("sphere1-mirror", lambda: _layered(_sphere, (0.0, 1.0), BC.reflective)),
    ("sphere3-mirror", lambda: _layered(_sphere, _MR3, BC.reflective)),
    ("sphere1-white", lambda: _layered(_sphere, (0.0, 1.0), BC.white)),
    ("slab3-mirrors", lambda: _layered(_slab, _SLB3, BC.reflective, BC.reflective)),
]
_CLOSED_MIXTURES = {"PU2": _PU2, "URRb": _URRB, "UP2N": _UP2N}
#: The exact infinite medium's gauge, written here: a fission production density of 100 n/cm^3/s
#: (``orpheus.derivations.common.exact_homogeneous``, the homogeneous solver's), not a finite body's.
_EXACT_DENSITY = 100.0


@pytest.mark.verifies("characteristic-door-gauge")
@pytest.mark.l1
@pytest.mark.verifies("characteristic-pencil")
@pytest.mark.parametrize("body", [b[1] for b in _CLOSED_BODIES], ids=[b[0] for b in _CLOSED_BODIES])
@pytest.mark.parametrize("name", list(_CLOSED_MIXTURES))
@pytest.mark.rests_on(_D5, _HERE + "test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux")
def test_a_closed_homogeneous_body_reads_the_exact_infinite_mediums_k_and_flux(body, name: str) -> None:
    """[K2; l1, ``characteristic-pencil``] A closed body of one mixture (under one material id per interval) and
    the exact infinite medium of that mixture (``exact_infinite_medium_reference``, rational arithmetic, a
    structurally independent reference): the same k to 1e-12 relative, and each group's flux integral of the
    weight 1 (one ``Symbolic`` observable, read by both references) in the ratio of the two gauges, to 1e-12.

    The body's flux is gauged by its total production, 1 over the body; the
    medium's by a production density of 100 per unit volume
    (:data:`_EXACT_DENSITY`). The flat body flux is then the medium's divided
    by 100 V, and its integral over the body is the medium's reading over 100.
    `[M]` 2026-10-07 (``ta/measure_5.log``): k <= 4.6e-14, flux <= 3.4e-14, both on
    URRb (upscatter). First reds: the gauge per unit
    volume (off by V); the gauge omitted (the unit-norm vector); the gauge's
    weight chi instead of nu Sigma_f.
    """
    mixture = _CLOSED_MIXTURES[name]
    geometry = body()
    spec = _spec(geometry, _same(mixture, len(geometry.mat_ids)), Eigen(_K))
    exact = exact_infinite_medium_reference(InfiniteMediumSpecification(0, mixture, Eigen(_K)))
    k = _read(spec, Eigenvalue())
    assert abs(k / exact.read(Eigenvalue()).value - 1.0) < 1e-12
    for g in range(2):
        weight = FluxIntegral(Symbolic.of(*(1 if h == g else 0 for h in range(2))))
        assert abs(_EXACT_DENSITY * _read(spec, weight) / exact.read(weight).value - 1.0) < 1e-12


def _productions(mat_ids, materials) -> tuple[np.ndarray, np.ndarray]:
    """Per region and group, in interval order: nu Sigma_f, and the default production nu Sigma_f + 2 sum_g' Sigma_2,g->g'
    (the (n,2n) multiplicity written here), from the mixtures' arrays."""
    fission = np.array([np.asarray(materials[m].SigP, dtype=float) for m in mat_ids])
    n2n = np.array([2.0 * materials[m].Sig2[0].toarray().sum(axis=1) for m in mat_ids])
    return fission, fission + n2n


_FISSION_GAUGE = CellCoefficient.every(Channel.FISSION_EMISSION)


@pytest.mark.verifies("characteristic-door-gauge")
@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux")
def test_a_heterogeneous_eigen_flux_produces_one_neutron() -> None:
    """[K3; l1] On the heterogeneous vacuum sphere (two fissile mixtures, one with an (n,2n) transfer, and an
    absorber; ids out of order), ``FluxIntegral`` of the default production, nu Sigma_f + 2 sum_g' Sigma_2
    per region and group written here from the mixtures, reads 1 (the declared default gauge, the user's
    ruling of 2026-10-08) to 1e-13; nu Sigma_f alone reads below 1 by more than 1e-4 (the activation: the
    (n,2n) emission is counted); every group's flux integral is positive.

    `[M]` 2026-10-08 (``ta/measure_8.log``): 0; nu Sigma_f alone reads 0.991. First reds: the gauge read from
    chi (the spectrum for the production); the gauge spelled fission alone
    whatever the question declares; the density gauge 100; the regions'
    production read in id order.
    """
    spec = _spec(_sphere(_MR3, (1, 0, 2), BC.vacuum), _WIRING_MATERIALS, Eigen(_K))
    fission, production = _productions((1, 0, 2), _WIRING_MATERIALS)
    assert abs(_read(spec, FluxIntegral(_rwc(production))) - 1.0) < 1e-13
    assert _read(spec, FluxIntegral(_rwc(fission))) < 1.0 - 1e-4
    for g in range(2):
        assert _read(spec, FluxIntegral(Symbolic.of(*(1 if h == g else 0 for h in range(2))))) > 0.0


@pytest.mark.verifies("characteristic-door-gauge")
@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_a_heterogeneous_eigen_flux_produces_one_neutron")
def test_a_declared_fission_gauge_rescales_the_flux_by_the_ratio_of_the_productions() -> None:
    """[K5; l1] The same body under the declared fission gauge: k is the default's bit for bit; nu Sigma_f
    reads 1 to 1e-13; and every flux integral (the two group totals, a region indicator) is the default's
    divided by the default answer's fission production, to 1e-13 (the declaration is read, and it changes the
    scale alone).

    `[M]` 2026-10-08 (``ta/measure_8.log``): k bitwise; 0; <= 1.1e-16. First reds: the declared gauge
    ignored (the door reads the default production whatever the question
    says, so the fission leg reads above 1); the gauge applied to k.
    """
    geometry = _sphere(_MR3, (1, 0, 2), BC.vacuum)
    default = _spec(geometry, _WIRING_MATERIALS, Eigen(_K))
    declared = _spec(geometry, _WIRING_MATERIALS, Eigen(_K, gauge=_FISSION_GAUGE))
    fission, _ = _productions((1, 0, 2), _WIRING_MATERIALS)
    assert _read(declared, Eigenvalue()) == _read(default, Eigenvalue())
    assert abs(_read(declared, FluxIntegral(_rwc(fission))) - 1.0) < 1e-13
    ratio = _read(default, FluxIntegral(_rwc(fission)))
    for weight in (Symbolic.of(1, 0), Symbolic.of(0, 1), _indicator(3, 1, 1).weight):
        assert abs(_read(declared, FluxIntegral(weight)) * ratio / _read(default, FluxIntegral(weight)) - 1.0) < 1e-13


@pytest.mark.l1
@pytest.mark.verifies("characteristic-pencil")
@pytest.mark.rests_on(_D2, _HERE + "test_the_doors_k_is_the_fundamental_of_the_system_written_by_hand")
def test_the_door_reads_k_one_at_soods_critical_sphere() -> None:
    """[K4; l1, ``characteristic-pencil``] Sood 2003 PU-2-0-SP (Siewert-Thomas 2G F_N): the vacuum sphere at its
    published critical radius, posed through a ``BC.vacuum`` tag, reads k = 1 within D2's band at the working
    point, 3.6e-6 (|dk/dmfp| x half a unit of the printed radius's last digit, measured in D2).

    `[M]` 2026-10-07: |k - 1| = 7.0e-7, D2's reading. First red: the vacuum tag
    read as a mirror (k moves to k_inf).
    """
    case = sood.PU_2_0_SP_STUB
    mixture = case.materials[0]
    published = case.truth.critical_dimension_mfp
    assert published is not None
    radius = published / float(mixture.SigT[0])
    spec = _spec(_sphere((0.0, radius), (0,), BC.vacuum), {0: mixture}, Eigen(_K))
    assert abs(_read(spec, Eigenvalue()) - 1.0) < 3.6e-6


# ── 5a.5 the mode nearest tau ────────────────────────────────────────────


def _k_of_buckling(mixture, buckling: float) -> float:
    r"""The k of the flux e^{iBx} in the infinite medium: (diag(Sigma_t / L(B)) - S)^-1 F's dominant eigenvalue,
    L(B) = arctan(B / Sigma_t) / (B / Sigma_t), the transport kernel's Fourier transform (isotropic emission)."""
    sigma_t = np.asarray(mixture.SigT, dtype=float)
    beta = buckling / sigma_t
    kernel = np.arctan(beta) / beta if buckling > 0.0 else np.ones_like(beta)
    matrix = np.linalg.solve(np.diag(sigma_t / kernel) - _scattering_by_hand(mixture), _fission_by_hand(mixture))
    return float(np.max(np.linalg.eigvals(matrix).real))


_WIDTH = 2.0


def _mirror_slab(mixture, mode: Any) -> GeometrySpecification:
    return _spec(_slab((0.0, _WIDTH), (0,), BC.reflective, BC.reflective), {0: mixture}, Eigen(_K, mode=mode))


@pytest.mark.l1
@pytest.mark.verifies("characteristic-pencil")
@pytest.mark.parametrize("name", ["PU2", "URRb"])
@pytest.mark.rests_on(_HERE + "test_the_doors_k_is_the_fundamental_of_the_system_written_by_hand")
def test_the_mode_nearest_tau_is_the_mirror_slabs_first_spatial_mode(name: str) -> None:
    """[N1; l1, ``characteristic-pencil``] A homogeneous slab between two mirrors (width 2 cm) has the exact
    eigenfunctions cos(n pi x / a): its mode n reads k(n pi / a), the dominant eigenvalue of
    (diag(Sigma_t / L(B)) - S)^-1 F with L(B) = arctan(B / Sigma_t) / (B / Sigma_t), a closed form written here.
    ``Nearest(tau)`` with tau 1 % above k(pi / a) reads it to 1e-6 relative, two orders apart from k_inf.

    A higher mode is not flux-shape independent (it is the transport
    kernel's response at B = pi / a), so this is the 2-group row the
    fundamental cannot be. `[M]` 2026-10-07: 6.5e-8 (PU2), 5.0e-7 (URRb) at the
    working point; the mode n = 2 is still 1e-3 off at degree 3 and is not
    asserted. First reds: ``Nearest`` read as ``Fundamental``; the spectrum
    searched by |k| order (the second element) on a body whose second element
    is not this mode (none here: the row is blind to that arm, declared).
    """
    mixture = _CLOSED_MIXTURES[name]
    k1 = _k_of_buckling(mixture, np.pi / _WIDTH)
    assert abs(k1 / _k_of_buckling(mixture, 0.0) - 1.0) > 0.5            # the activation: the mode is not k_inf
    k = _read(_mirror_slab(mixture, Nearest(1.01 * k1)), Eigenvalue())
    assert abs(k / k1 - 1.0) < 1e-6


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_mode_nearest_tau_is_the_mirror_slabs_first_spatial_mode")
def test_tau_is_read_in_the_k_chart() -> None:
    """[N2; l1] With tau = 1, between the harmonic and the arithmetic mean of k_inf (2.68) and k(pi / a) (0.286)
    for PU2, the nearest mode in k is k(pi / a); in the 1/k chart it would be k_inf.

    `[M]` 2026-10-07: reads k(pi / a). First red: the distance taken as
    |1/k - 1/tau|.
    """
    k0, k1 = _k_of_buckling(_PU2, 0.0), _k_of_buckling(_PU2, np.pi / _WIDTH)
    tau = 1.0
    assert 2.0 * k0 * k1 / (k0 + k1) < tau < 0.5 * (k0 + k1)               # the premise: the charts disagree at tau
    assert abs(_read(_mirror_slab(_PU2, Nearest(tau)), Eigenvalue()) / k1 - 1.0) < 1e-6


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_mode_nearest_tau_is_the_mirror_slabs_first_spatial_mode")
def test_tau_near_the_fundamental_reads_the_fundamentals_k() -> None:
    """[N3; l1] ``Nearest(0.99 k_inf)`` reads the fundamental's k to 1e-13 relative on the heterogeneous vacuum
    sphere. The two eigenvalues come from two eigensolver calls (``fundamental`` and ``spectrum``).
    `[M]` 2026-10-07: 0. First red: the Nearest search over the spectrum's absolute values in the 1/k chart
    reaching another mode (none here: blind, declared; N2 carries the chart)."""
    fundamental = _spec(_sphere(_MR3, (1, 0, 2), BC.vacuum), _WIRING_MATERIALS, Eigen(_K))
    k = _read(fundamental, Eigenvalue())
    nearest = _spec(_sphere(_MR3, (1, 0, 2), BC.vacuum), _WIRING_MATERIALS, Eigen(_K, mode=Nearest(0.99 * k)))
    assert abs(_read(nearest, Eigenvalue()) / k - 1.0) < 1e-13


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_mode_nearest_tau_is_the_mirror_slabs_first_spatial_mode")
def test_tau_zero_reads_a_genuine_mode_not_the_null_fission_cluster() -> None:
    """[N5; l1] ``Nearest(0)`` on the mirror slab (PU2: chi (x) nu Sigma_f has rank one per node, so the pencil
    carries a cluster of eigenvalues at k = 0 whose vectors fission does not see) reads a real eigenvalue above
    1e-6 and below k(pi / a): the smallest genuine mode, never the cluster's rounding residue.

    One-sided, as the claim is: the smallest resolved spatial mode moves
    with the resolution, the cluster sits at rounding. `[M]` 2026-10-07
    (``ta/measure_7.log``): 4.3e-3 at the working point; the cluster's
    eigenvalues are 2.8e-18 and up. First red: the null-production filter removed
    (the nearest eigenvalue is then a residue of order 1e-16, or complex and
    refused).
    """
    k = _read(_mirror_slab(_PU2, Nearest(0.0)), Eigenvalue())
    assert 1e-6 < k < _k_of_buckling(_PU2, np.pi / _WIDTH), k


_HIGHER_MODE_REFUSAL = "a higher mode has no flux scale"


@pytest.mark.foundation
@pytest.mark.parametrize("at", ["the-first-spatial-mode", "the-fundamental"])
@pytest.mark.rests_on(_HERE + "test_the_mode_nearest_tau_is_the_mirror_slabs_first_spatial_mode")
def test_a_nearest_answer_reads_its_eigenvalue_and_refuses_a_flux_integral(at: str) -> None:
    """[N4] A ``Nearest`` answer reads ``Eigenvalue`` only: a ``FluxIntegral`` of it, and a ``Ratio`` of two,
    are refused with NotImplementedError ("a higher mode has no flux scale", a SCOPE-BOUNDARY: the user's
    ruling of 2026-10-07), at tau on the mirror slab's first spatial mode AND at tau near the fundamental
    (the refusal is on the question's mode, not on the value it finds).

    Why (`[M]` 2026-10-07, ``ta/probe_n4.py``): every higher mode of a
    closed homogeneous body is biorthogonal to the flat adjoint, so its net
    fission production is zero, and the production gauge (then 100) divided by a
    rounding residue (-1.3e-16 at degree 3): the left-half group-0 integral
    of the mirror slab's first mode read -8.4e16 (-1.0e17 at degree 4) before
    the ruling. First reds: the refusal keyed on the mode found rather than
    on the question (the fundamental leg constructs a reading); the refusal
    removed (the first leg reads a number).
    """
    k1 = _k_of_buckling(_PU2, np.pi / _WIDTH)
    tau = k1 if at == "the-first-spatial-mode" else 0.99 * _k_of_buckling(_PU2, 0.0)
    spec = _mirror_slab(_PU2, Nearest(tau))
    assert _read(spec, Eigenvalue()) > 0.0
    left_half = FluxIntegral(Symbolic.of(sympy.Piecewise((1, Symbolic.r < sympy.Rational(1)), (0, True)), 0))
    with pytest.raises(NotImplementedError, match=_HIGHER_MODE_REFUSAL):
        _read(spec, left_half)
    with bypass(), pytest.raises(NotImplementedError, match=_HIGHER_MODE_REFUSAL):
        characteristic_reference(spec, _RES).read(Ratio(left_half, FluxIntegral(Symbolic.of(1, 0))))


# ── 5a.6 the fixed source ────────────────────────────────────────────────


_E2_BODIES = [
    pytest.param("sphere", _MR3, lambda: _layered(_sphere, _MR3, BC.reflective), id="sphere3-mirror"),
    pytest.param("slab", _SLB3, lambda: _layered(_slab, _SLB3, BC.reflective, BC.reflective), id="slab3-mirrors"),
    pytest.param("sphere", (0.4, 2.0), lambda: _sphere((0.4, 2.0), (0,), BC.reflective, BC.white), id="hollow-white-mirror"),
    pytest.param("cylinder", (0.0, 1.0), lambda: StructuredGeometry.cylinder((0.0, 1.0), (0,), outer=BC.reflective),
                 marks=pytest.mark.slow, id="cylinder1-mirror"),
]


@pytest.mark.l1
@pytest.mark.verifies("characteristic-fixed-source")
@pytest.mark.parametrize(("chart", "breakpoints", "body"), _E2_BODIES)
@pytest.mark.rests_on(_E2, _HERE + "test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux")
def test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux_through_the_door(chart, breakpoints, body) -> None:
    """[F1; l1, ``characteristic-fixed-source``] E2 through the door: a closed body of ``_UP2N`` (upscatter,
    (n,2n), subcritical fission) with the region-wise source q = (1, 0.4) everywhere; each region's group flux
    integral divided by the region's volume (written here) is (diag Sigma_t - S - F)^-1 q, to 1e-13.

    The region-wise source is a rate, read as given. `[M]` 2026-10-07: <= 5.8e-15;
    the cylinder 2.2e-14 (slow, 36 s). First reds: the region-wise
    source read with the retraction's 4 pi; the indicator weights read in id
    order on a body whose ids are the interval order (blind, declared: the
    mixture is one).
    """
    geometry = body()
    n = len(geometry.mat_ids)
    q = np.array([1.0, 0.4])
    spec = _spec(geometry, _same(_UP2N, n), FixedSource(_rwc(np.tile(q, (n, 1)))))
    loss, production = _zero_d(_UP2N)
    expected = np.linalg.solve(loss - production, q)
    volumes = _volumes(chart, breakpoints)
    for region in range(n):
        for g in range(2):
            reading = _read(spec, _indicator(n, region, g)) / volumes[region]
            assert abs(reading / expected[g] - 1.0) < 1e-13, (region, g, reading, expected[g])


@pytest.mark.l1
@pytest.mark.rests_on(_SU2, _HERE + "test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux")
def test_a_source_where_its_group_emits_nothing_widens_the_support_and_balances() -> None:
    """[F2; l1] A closed mirror sphere (``_ABS`` | ``_SELF_ONLY`` | ``_ABS``) with a group-1 source only in the
    middle region, where group 1 is void and nothing is emitted into it: the door poses the source's region, and
    the absorption read through ``FluxIntegral`` (Sigma_a written here) equals the source's rate times the
    region's volume, to 1e-12.

    A balance is blind to a redistribution that conserves (``vv``
    anti-#8): this row carries the support's widening, the value rows carry
    the rest. `[M]` 2026-10-07: 1.1e-15. First red: the system posed without
    the source's regions (``EmissionSpace.restrict`` refuses the source).
    """
    mixtures = {0: _ABS, 1: _SELF_ONLY}
    source = np.array([[0.0, 0.0], [0.0, 1.5], [0.0, 0.0]])
    spec = _spec(_sphere(_MR3, (0, 1, 0), BC.reflective), mixtures, FixedSource(_rwc(source)))
    absorption = np.array([np.asarray(mixtures[m].SigT) - mixtures[m].SigS[0].toarray().sum(axis=1) for m in (0, 1, 0)])
    rate = 1.5 * _volumes("sphere", _MR3)[1]
    assert abs(_read(spec, FluxIntegral(_rwc(absorption))) / rate - 1.0) < 1e-12


@pytest.mark.foundation
@pytest.mark.parametrize("question", [FixedSource(_rwc(np.ones((3, 2)))), Response(_rwc(np.ones((3, 2))))],
                         ids=["source", "response"])
@pytest.mark.rests_on(_E6)
def test_a_source_question_on_a_body_supercritical_by_n2n_is_refused_through_the_door(question) -> None:
    """[F3] E6 through the door: the closed mirror sphere of ``_N2N`` (supercritical by (n,2n) alone) refuses a
    source and a detector with ``NoLeastSolution`` ("spectral radius of the gain") at the first reading.
    The transposed problem has the same gain spectrum, so the response is refused too. First red: the source
    question answered on the k pencil (whose gain is the fission alone, none here)."""
    spec = _spec(_layered(_sphere, _MR3, BC.reflective), _same(_N2N, 3), question)
    with pytest.raises(NoLeastSolution, match="spectral radius of the gain"):
        _read(spec, _indicator(3, 0, 0))


_OPEN_HETERO = _sphere(_ADJ_BREAKPOINTS, (0, 1, 2), BC.vacuum)
_ADJ_MATERIALS = {0: _ADJ_FUEL, 1: _ADJ_ABSORBER, 2: _ADJ_REFLECTOR}


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux",
                      _HERE + "test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux_through_the_door")
@pytest.mark.verifies("characteristic-door-lifts")
def test_a_symbolic_source_is_a_density_over_directions_and_a_symbolic_weight_or_detector_is_not() -> None:
    """[F4; l1] On the heterogeneous vacuum sphere (qa's adjoint body), against the region-wise readings:
    (a) a constant ``Symbolic`` source q reads as the region-wise source 4 pi q; (b) a constant ``Symbolic``
    weight w reads as the region-wise weight w; (c) a constant ``Symbolic`` detector d reads as the region-wise
    detector d (both retract to 4 pi d). Each to 1e-13 relative.

    The 4 pi is the retraction of a density over the sphere of directions
    (``angular_measure.retraction``), not a typed constant. `[M]` 2026-10-07:
    2.2e-16, 0 and 2.2e-16. First reds, one per leg: the source read
    by the pullback (a, off by 4 pi); the weight read by the retraction (b);
    the symbolic detector retracted and the region-wise one not (c).
    """
    q, w, d = np.array([0.8, 0.3]), np.array([0.5, 1.7]), np.array([1.2, 0.4])
    uniform = lambda v: _rwc(np.tile(v, (3, 1)))                          # noqa: E731
    symbolic_source = _read(_spec(_OPEN_HETERO, _ADJ_MATERIALS, FixedSource(_symbolic_constant(q))), FluxIntegral(uniform(w)))
    regionwise_source = _read(_spec(_OPEN_HETERO, _ADJ_MATERIALS, FixedSource(uniform(4.0 * np.pi * q))), FluxIntegral(uniform(w)))
    assert abs(symbolic_source / regionwise_source - 1.0) < 1e-13
    source = _spec(_OPEN_HETERO, _ADJ_MATERIALS, FixedSource(uniform(q)))
    assert abs(_read(source, FluxIntegral(_symbolic_constant(w))) / _read(source, FluxIntegral(uniform(w))) - 1.0) < 1e-13
    symbolic_detector = _read(_spec(_OPEN_HETERO, _ADJ_MATERIALS, Response(_symbolic_constant(d))), FluxIntegral(uniform(q)))
    regionwise_detector = _read(_spec(_OPEN_HETERO, _ADJ_MATERIALS, Response(uniform(d))), FluxIntegral(uniform(q)))
    assert abs(symbolic_detector / regionwise_detector - 1.0) < 1e-13


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux_through_the_door")
def test_a_regionwise_source_is_read_onto_the_nodes_whatever_the_projection_points() -> None:
    """[F5; l1] A region-wise source and weight are read onto the nodes exactly (the basis is nodal, the
    fix (1) of 2026-10-07: ``PanelBasis.on_nodes``), so the readings at ``source_points`` 8 and 10
    are equal bit for bit, and a source's regions are its non-zero rows: a source zero in the middle region
    poses no source there (the system's ``source_regions`` row is False wherever the table's row is zero).

    `[M]` 2026-10-07: equal bitwise after fix (1); before it the two readings
    differed in the last bits (the L2 projection of a constant is exact to
    rounding, not bitwise). First reds: the
    region-wise table projected through the quadrature; the source regions
    read from the projected coefficients.
    """
    table = np.array([[1.0, 0.4], [0.0, 0.0], [0.7, 0.2]])
    reads = []
    for points in (8, 10):
        resolution = _resolution_with(source_points=points)
        spec = _spec(_OPEN_HETERO, _ADJ_MATERIALS, FixedSource(_rwc(table)))
        reads.append(_read(spec, FluxIntegral(_rwc(np.full((3, 2), 0.5))), resolution))
        regions = _derivation(spec, resolution).system.source_regions
        assert regions is not None and not regions[1].any()
    assert reads[0] == reads[1], reads


# ── 5a.7 the response ────────────────────────────────────────────────────


_RECIPROCITY_BODIES = [
    ("sphere-vacuum", lambda: _sphere(_ADJ_BREAKPOINTS, (0, 1, 2), BC.vacuum)),
    ("slab-mirror-vacuum", lambda: _slab(_ADJ_BREAKPOINTS, (0, 1, 2), BC.reflective, BC.vacuum)),
    ("sphere-white", lambda: _sphere(_ADJ_BREAKPOINTS, (0, 1, 2), BC.white)),
]


@pytest.mark.verifies("characteristic-door-reciprocity")
@pytest.mark.l1
@pytest.mark.verifies("characteristic-adjoint")
@pytest.mark.parametrize("body", [b[1] for b in _RECIPROCITY_BODIES], ids=[b[0] for b in _RECIPROCITY_BODIES])
@pytest.mark.rests_on(_D13IV, _HERE + "test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux")
def test_a_detectors_reading_of_a_source_is_the_sources_reading_of_the_detectors_importance(body) -> None:
    """[A1; l1, ``characteristic-adjoint``] Reciprocity across two questions: ``FixedSource(q)`` read by
    ``FluxIntegral(r)`` equals ``Response(r)`` read by ``FluxIntegral(q)`` over 4 pi, two seeded positive
    region-wise pairs, to 1e-13 relative, on qa's adjoint body (up- and downscatter, (n,2n), fission with chi in both
    groups, supports that differ by group).

    The response's answer is the adjoint SCALAR flux R psi^dagger (the
    ruling of 2026-10-07): the detector is pulled back by R^dagger and the
    transposed problem's source is its retraction, 4 pi r for a region-wise
    r (the 4 pi written here, the solid angle). The two sides are two systems: the forward on the cross sections with
    q's regions, the response on the TRANSPOSED cross sections with r's
    regions; their emission spaces differ. Declared stabiliser: both share
    the transport blocks K_g (one Sigma_t), so an error in K symmetric under
    i <-> j is invisible here (D13 (iv) and the rung-3 rows carry K). The
    reading's premise is asserted: the scattering and fission tables are
    not symmetric. `[M]` 2026-10-07 (``ta/measure_5.log``): <= 7.8e-16. First reds:
    ``transposed()`` returning the tables untransposed; transposing the
    scattering only; the detector pulled back and not retracted (off by 4 pi).
    """
    geometry = body()
    xs = RegionCrossSections.of([_ADJ_FUEL, _ADJ_ABSORBER, _ADJ_REFLECTOR])
    assert not np.array_equal(xs.scattering[0], xs.scattering[0].T) and not np.array_equal(xs.fission[0], xs.fission[0].T)
    rng = np.random.default_rng(7)
    for _ in range(2):
        q, r = _rwc(rng.random((3, 2))), _rwc(rng.random((3, 2)))
        forward = _read(_spec(geometry, _ADJ_MATERIALS, FixedSource(q)), FluxIntegral(r))
        adjoint = _read(_spec(geometry, _ADJ_MATERIALS, Response(r)), FluxIntegral(q))
        assert abs(adjoint / (4.0 * np.pi * forward) - 1.0) < 1e-13, (forward, adjoint)


@pytest.mark.verifies("characteristic-door-response")
@pytest.mark.l1
@pytest.mark.verifies("characteristic-adjoint")
@pytest.mark.parametrize(("chart", "breakpoints", "body"), _E2_BODIES[:2])
@pytest.mark.rests_on(_HERE + "test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux_through_the_door")
def test_a_closed_bodys_importance_is_the_infinite_mediums_adjoint_solve(chart, breakpoints, body) -> None:
    """[A2; l1, ``characteristic-adjoint``] A closed body of ``_UP2N`` and the uniform detector r = (0.3, 1):
    each region's group integral of the adjoint scalar flux, over the region's volume, is
    4 pi (diag Sigma_t - S - F)^-T r (the detector's retraction as the source), to 1e-13; the untransposed solve differs from it by more than 1e-2 (the premise: the transpose matters).

    `[M]` 2026-10-07 (``ta/measure_5.log``): <= 6.2e-15; the untransposed solve 1.5 off.
    First reds: ``transposed()`` returning the tables untransposed; the
    fission left untransposed.
    """
    geometry = body()
    n = len(geometry.mat_ids)
    r = np.array([0.3, 1.0])
    spec = _spec(geometry, _same(_UP2N, n), Response(_rwc(np.tile(r, (n, 1)))))
    loss, production = _zero_d(_UP2N)
    expected = 4.0 * np.pi * np.linalg.solve((loss - production).T, r)
    untransposed = 4.0 * np.pi * np.linalg.solve(loss - production, r)
    assert np.max(np.abs(untransposed / expected - 1.0)) > 1e-2
    volumes = _volumes(chart, breakpoints)
    for region in range(n):
        for g in range(2):
            reading = _read(spec, _indicator(n, region, g)) / volumes[region]
            assert abs(reading / expected[g] - 1.0) < 1e-13, (region, g, reading, expected[g])


@pytest.mark.l1
@pytest.mark.verifies("characteristic-adjoint")
@pytest.mark.rests_on(_SYSTEM + "test_a_detectors_reading_of_a_source_is_its_response_paired_with_the_source",
                      _HERE + "test_a_detectors_reading_of_a_source_is_the_sources_reading_of_the_detectors_importance")
@pytest.mark.verifies("characteristic-door-response")
def test_the_doors_response_is_the_forward_systems_adjoint_on_the_same_posing() -> None:
    """[A3; l1, ``characteristic-adjoint``] The door's ``Response(r)`` read by ``FluxIntegral(q)`` equals
    4 pi ``system.response(r)`` paired with W_s q (``system.response`` is the importance per unit emission rate,
    the door's answer the adjoint scalar flux) on the FORWARD system posed here by hand (mixtures in interval order,
    walls as ``Wall`` tuples, q's regions as its source regions), to 1e-13 relative: the door's route (the
    transposed problem's forward solve) against rung 4's (the adjoint pencil of the forward system).

    The two share the transport blocks K_g and nothing above them: one
    solves (W_s' - (S' + F')K') q' = W_s' r on the transposed supports, the
    other (W_s - K^T (S + F)^T) phi^dagger = K^T r on the forward ones.
    `[M]` 2026-10-07: 0 (``ta/measure_5.log``). First reds: ``transposed()``
    returning the tables untransposed; the detector not retracted (off by 4 pi).
    """
    q_table, r_table = np.array([[1.0, 0.2], [0.0, 0.5], [0.3, 0.0]]), np.array([[0.2, 1.0], [0.6, 0.0], [0.0, 0.9]])
    geometry = _sphere(_ADJ_BREAKPOINTS, (0, 1, 2), BC.vacuum)
    door = _read(_spec(geometry, _ADJ_MATERIALS, Response(_rwc(r_table))), FluxIntegral(_rwc(q_table)))
    basis = _basis("sphere", _ADJ_BREAKPOINTS)
    system = GalerkinSystem(basis, _walls("sphere", _ADJ_BREAKPOINTS, (_VACUUM,)),
                            RegionCrossSections.of([_ADJ_FUEL, _ADJ_ABSORBER, _ADJ_REFLECTOR]), TransportResolution(8, 12, 12),
                            source_regions=q_table != 0.0)
    q, r = q_table[basis.region].T, r_table[basis.region].T
    by_hand = float(system.response(r) @ system.emission_mass @ system.emission.restrict(q))
    assert abs(door / (4.0 * np.pi * by_hand) - 1.0) < 1e-13, (door, by_hand)


# ── 5a.8 the reference solution, content identity and the traced memo ────


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux")
def test_the_factory_returns_an_uncertified_reference_that_reads_a_ratio_as_a_quotient() -> None:
    """[R3] ``characteristic_reference`` returns a ``ReferenceSolution`` with no certificate around a
    ``CharacteristicDerivation`` of the same specification and resolution; ``read(Ratio(a, b))`` is the
    quotient of the two readings, bit for bit; an eigen ratio is independent of the gauge (read against
    the unscaled Galerkin flux here, 1e-13). First red: the eigen flux gauged per observable."""
    spec = _spec(_sphere(_MR3, (1, 0, 2), BC.vacuum), _WIRING_MATERIALS, Eigen(_K))
    solution = characteristic_reference(spec, _RES)
    assert type(solution) is ReferenceSolution and solution.certificate is None
    assert isinstance(solution.derivation, CharacteristicDerivation)
    assert solution.derivation == CharacteristicDerivation(spec, _RES)
    a, b = _indicator(3, 0, 0), _indicator(3, 2, 1)
    with bypass():
        quotient = solution.read(Ratio(a, b)).value
        assert quotient == solution.read(a).value / solution.read(b).value
    derivation = _derivation(spec)
    flux = derivation.system.flux(derivation.system.pencil.fundamental().vector)
    by_hand = _volume_integral(derivation.basis, flux, lambda c: np.stack([c < _MR3[1], 0 * c]).astype(float), None) \
        / _volume_integral(derivation.basis, flux, lambda c: np.stack([0 * c, c > _MR3[2]]).astype(float), None)
    assert abs(quotient / by_hand - 1.0) < 1e-13


def _derivation_with(**changes: Any) -> CharacteristicDerivation:
    fields = {"specification": _spec(_sphere(_MR3, (1, 0, 2), BC.vacuum), _WIRING_MATERIALS, Eigen(_K)),
              "resolution": _RES} | changes
    return CharacteristicDerivation(**fields)


def _resolution_with(**changes: Any) -> Resolution:
    """The working point with 10 projection points, enough for degree 4 (so the degree leg moves one field)."""
    return Resolution(**({"degree": 3, "layers": 2, "ratio": 0.4, "transport": TransportResolution(8, 12, 12),
                          "source_points": 10} | changes))


ROSTER: tuple[Entry, ...] = (
    Entry(
        cls=TransportResolution,
        base=lambda: TransportResolution(8, 12, 12),
        parts=("line_points", "points", "inner_points"),
        perturb={
            "line_points": (leg("more line points", lambda: TransportResolution(9, 12, 12)),),
            "points": (leg("more traversal points", lambda: TransportResolution(8, 13, 12)),),
            "inner_points": (leg("more inner points", lambda: TransportResolution(8, 12, 13)),),
        },
    ),
    Entry(
        cls=Resolution,
        base=_resolution_with,
        parts=("degree", "layers", "ratio", "transport", "source_points"),
        perturb={
            "degree": (leg("another degree", lambda: _resolution_with(degree=4)),),
            "layers": (leg("another grading depth", lambda: _resolution_with(layers=3)),),
            "ratio": (leg("another grading ratio", lambda: _resolution_with(ratio=0.3)),),
            "transport": (leg("another transport", lambda: _resolution_with(transport=TransportResolution(8, 12, 13))),),
            "source_points": (leg("more projection points", lambda: _resolution_with(source_points=11)),),
        },
    ),
    Entry(
        cls=CharacteristicDerivation,
        base=_derivation_with,
        parts=("specification", "resolution"),
        perturb={
            "specification": (leg("another body", lambda: _derivation_with(
                specification=_spec(_sphere((0.0, 0.5, 1.5, 2.1), (1, 0, 2), BC.vacuum), _WIRING_MATERIALS, Eigen(_K)))),),
            "resolution": (leg("another resolution", lambda: _derivation_with(resolution=_resolution_with(degree=4))),),
        },
        pairs=(("a-fresh-specification-and-resolution", _derivation_with,
                lambda: _derivation_with(specification=_spec(_sphere(_MR3, (1, 0, 2), BC.vacuum), dict(_WIRING_MATERIALS),
                                                             Eigen(_K)), resolution=Resolution(3, 2, 0.4, TransportResolution(8, 12, 12), 8))),),
    ),
)


@pytest.mark.foundation
@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_c1_the_population_is_the_parts(entry: Entry) -> None:
    """[C1] Each value's content is its constructor fields; the derivation's system, basis, cross sections
    and source are derived and not content. First red: a derived field compared."""
    check_population(entry)


@pytest.mark.foundation
@pytest.mark.parametrize("entry, part, the_leg", perturbation_ids(ROSTER),
                         ids=[param_id(e.id, p, lg[0]) for e, p, lg in perturbation_ids(ROSTER)])
def test_c1_each_part_moves_the_digest(entry: Entry, part: str, the_leg: Any) -> None:
    """[C1] Every constructor field moves the digest. First red: a field excluded from the content."""
    check_perturbation(entry, part, the_leg)


@pytest.mark.foundation
@pytest.mark.parametrize("entry, pair", pair_ids(ROSTER), ids=[param_id(e.id, p[0]) for e, p in pair_ids(ROSTER)])
def test_c1_equal_fields_are_one_value(entry: Entry, pair: Any) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.foundation
def test_c2_evaluate_is_a_traced_memo() -> None:
    """[C2] ``CharacteristicDerivation.evaluate`` is bound as a ``TracedMemo``. First red: the decorator dropped."""
    assert isinstance(CharacteristicDerivation.__dict__["evaluate"], memo_api.name("TracedMemo"))


def _tiny_spec(question: Any = None) -> GeometrySpecification:
    return _spec(_sphere((0.0, 1.0), (0,), BC.vacuum), {0: _UP2N}, Eigen(_K) if question is None else question)


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_c2_evaluate_is_a_traced_memo",
                      "tests/gates/numerics/test_traced_memo_process.py::test_m3_1_a_miss_generates_once_and_a_hit_never")
def test_c3_an_equal_derivation_reads_from_the_memo_and_a_new_resolution_does_not(tmp_path, monkeypatch) -> None:
    """[C3] A first reading of a one-region sphere at the smallest resolution starts an interpreter; a FRESH
    derivation of equal content reads the same observable with 0 interpreters started and the same value
    bit for bit; a derivation at another resolution is ``Absent`` (its key differs). First reds: the resolution
    left out of the content (the last leg reads ``Hit``); the evaluate undecorated (the first count is 0)."""
    spawns = memo_api.SpawnCounter(monkeypatch)
    with memo_api.cache_root(tmp_path):
        cold = CharacteristicDerivation(_tiny_spec(), _TINY).evaluate(Eigenvalue())
        started = spawns.count
        warm = CharacteristicDerivation(_tiny_spec(), _TINY).evaluate(Eigenvalue())
        other = CharacteristicDerivation(_tiny_spec(), _resolution_with(degree=2, layers=0, transport=TransportResolution(4, 4, 4),
                                                                         source_points=6))
        verdict = CharacteristicDerivation.__dict__["evaluate"].lookup(other, Eigenvalue())
    assert started >= 1, "the activation leg: the cold read generated nothing"
    assert spawns.count == started, f"the warm read started {spawns.count - started} interpreters"
    assert float(warm.value).hex() == float(cold.value).hex()
    assert memo_api.verdict_kind(verdict) == "Absent"


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_c3_an_equal_derivation_reads_from_the_memo_and_a_new_resolution_does_not")
def test_c4_a_refusal_crosses_the_memo_with_its_type(tmp_path) -> None:
    """[C4] A source on a body supercritical by (n,2n), read through the memo (a generating interpreter),
    raises ``NoLeastSolution`` here with its message. First red: the child's failure re-raised as another type."""
    spec = _spec(_sphere((0.0, 1.0), (0,), BC.reflective), {0: _N2N}, FixedSource(_rwc(np.ones((1, 2)))))
    with memo_api.cache_root(tmp_path), pytest.raises(NoLeastSolution, match="spectral radius of the gain"):
        CharacteristicDerivation(spec, _TINY).evaluate(_indicator(1, 0, 0))
