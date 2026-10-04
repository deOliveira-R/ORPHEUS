r"""The trajectory-resolvent reference as a ``ReferenceSolution`` (#405 P2 step 7b.2.2, gates R7b2.2–R7b2.7).

The multi-region trajectory resolvent (Variant α) becomes a reference
solution: lazy, with no certificate (its family derives no bound, #566 and,
for the cylinder, #516), every reading ``Uncertified``, read as the transport
integral of its emission density (the user's ruling 1 of 2026-10-03: one
reference, one answer), through ONE transport: the chord oracle with its
evaluation radii (``at=``) separate from the spline's knots. Spec §1.7b.2,
"The trajectory-resolvent references as ``ReferenceSolution``s" and "The
gates". Every production name goes through
:mod:`tests.gates.derivations._trajectory_resolvent_api`.

First red on ``65b0de93``: ``ModuleNotFoundError`` for the reference module
(R7b2.2–R7b2.6, R7b2.7's module rows); the oracle rows refuse the ``at=``
keyword (``TypeError``); ``Symbolic`` has no ``steps``.

Fixtures, at the gates' resolutions (the sphere's rows take seconds; the cylinder's 10–125 s, mostly the
reading's flux integrals, and are ``slow`` except one cheap reading law, R7b2.2.1):

* the A|B|A sphere and cylinder (:mod:`._aba_reference`, isotropic mixtures);
* the UNIFORM two-region sphere and cylinder: fuel A under two material ids,
  so the body is layered (the multi-region route) and its flux is flat in
  space and angle (the V_α1 identity): every point reads one value per group,
  the group ratio is the infinite medium's, exact (``exact_infinite_medium``,
  rational arithmetic, a structurally independent ground), and a cell
  indicator over the whole body reads its measure's fraction exactly.

Declared stabiliser: an eigen reference's readings are meaningful only as
ratios (its representative's scale is its own gauge), so a uniform scaling of
the emission density (a dropped 1/(4π)) is invisible to every gate here by
design, and harmless to every ratio a comparison reads.
"""
from __future__ import annotations

import ast
import functools
import importlib
import math
import pathlib
from fractions import Fraction
from typing import Any

import numpy as np
import pytest
import sympy

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.derivations.common.exact_homogeneous import exact_infinite_medium_of
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
from orpheus.numerics.observable import Eigenvalue, FluxIntegral, PointValue, Ratio
from orpheus.numerics.outcome import Measured, NotYet
from orpheus.numerics.question import Eigen, FixedSource
from orpheus.numerics.traced_memo import bypass
from orpheus.reference.reading import Uncertified
from orpheus.reference.verification import ReferenceNotValid, compare_uncertified, verify_agreement, verify_order
from orpheus.specification.specification import GeometrySpecification, InfiniteMediumSpecification
from tests.gates._content_identity_helpers import (
    Entry,
    check_equal_pair,
    check_perturbation,
    check_population,
    leg,
    pair_ids,
    param_id,
    perturbation_ids,
)
from tests.gates.derivations import _trajectory_resolvent_api as api
from tests.gates.sn.verification.analytical._aba_reference import aba_specification, isotropic_mixture

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/derivations/test_trajectory_resolvent_reference.py"
_FIRST_LEG = "tests/gates/derivations/test_trajectory_resolvent_regionwise_source.py::test_mr_oracle_first_leg_matches_the_line_integral"
_COORDS = (CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL)
_IDS = lambda c: c.name.lower()  # noqa: E731
_K = CellCoefficient.every(Channel.FISSION_EMISSION)
#: The coordinate systems for rows that read: the cylinder's rows cost 10–125 s at the gates' resolution (``[M]``
#: 2026-10-03, ``--durations=0``; the solve about 19 s, the rest the reading), so its case is ``slow``; the sphere
#: keeps every law in ``not slow``, and R7b2.2.1 keeps one cylinder reading law there at a cheaper resolution.
_COORDS_SOLVING = (CoordSystem.SPHERICAL, pytest.param(CoordSystem.CYLINDRICAL, marks=pytest.mark.slow, id="cylindrical"))


@functools.cache
def _uniform_specification(coord: CoordSystem) -> GeometrySpecification:
    """Fuel A under material ids 0 and 1 (regions [0, 1] and [1, 2] cm): a layered body with a flat solution."""
    fuel = isotropic_mixture("A")
    geometry = StructuredGeometry(coord=coord, breakpoints=(0.0, 1.0, 2.0), mat_ids=(0, 1), boundaries=(BC.reflective,))
    return GeometrySpecification(Materials({0: fuel, 1: fuel}), geometry, Eigen(_K))


@functools.cache
def _uniform(coord: CoordSystem):
    return api.reference(_uniform_specification(coord))


@functools.cache
def _aba(coord: CoordSystem):
    return api.reference(aba_specification(coord))


def _indicator(a: float, b: float, group: int) -> FluxIntegral:
    low, high = (sympy.Rational(str(x)) for x in (a, b))  # the decimal as written: 0.6 is 3/5
    step = sympy.Piecewise((1, (Symbolic.r >= low) & (Symbolic.r < high)), (0, True))
    return FluxIntegral(Symbolic.of(*(step if g == group else 0 for g in range(2))))


def _group_total(group: int) -> FluxIntegral:
    return FluxIntegral(Symbolic.of(*(1 if g == group else 0 for g in range(2))))


class _SolveSpy:
    """Counts calls to the two multi-region solvers, rebinding every module attribute that names them.

    A test using it reads :func:`in_process`, so the solves it counts run where it counts them.
    """

    def __init__(self, monkeypatch: pytest.MonkeyPatch) -> None:
        self.calls = 0
        for module_name, name in api.SOLVERS:
            module = importlib.import_module(module_name)
            original = getattr(module, name)

            def counting(*args, _original=original, **kwargs):
                self.calls += 1
                return _original(*args, **kwargs)

            monkeypatch.setattr(module, name, counting)


@pytest.fixture
def in_process():
    """Every memoised call runs in this process for the test (:func:`~orpheus.numerics.traced_memo.bypass`):
    a spy counts calls here, and a memoised reading would run its solve and its weight's steps in a generating
    process no spy reaches, where it counts 0 whether or not they ran (#405 P3)."""
    with bypass():
        yield


# ─────────────────────────────────────────────────────────────────────
# R7b2.2 — the factory: a reference with no certificate, solved lazily, refusals
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("coord", _COORDS, ids=_IDS)
def test_r7b2_2_the_factory_returns_an_uncertified_reference(coord) -> None:
    from orpheus.reference.solution import ReferenceSolution

    ref = api.reference(aba_specification(coord))
    assert type(ref) is ReferenceSolution
    assert ref.certificate is None
    assert ref.specification == aba_specification(coord)
    assert isinstance(ref.derivation, api.derivation_class())


@pytest.mark.parametrize("coord", _COORDS_SOLVING, ids=_IDS)
@pytest.mark.usefixtures("in_process")
def test_r7b2_2_construction_solves_nothing_and_the_solve_is_cached(coord, monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE: building the reference calls no solver; the first reading solves once; a second reading reuses it."""
    spy = _SolveSpy(monkeypatch)
    ref = api.reference(aba_specification(coord))
    assert spy.calls == 0, "the factory solved eagerly"
    ref.read(Eigenvalue())
    assert spy.calls == 1
    ref.read(Eigenvalue())
    ref.read(_group_total(0))
    assert spy.calls == 1, "the solve was not cached on the derivation"


def _slab_layered() -> GeometrySpecification:
    fuel = isotropic_mixture("A")
    geometry = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 1), left=BC.reflective, right=BC.reflective)
    return GeometrySpecification(Materials({0: fuel, 1: isotropic_mixture("B")}), geometry, Eigen(_K))


def _anisotropic() -> GeometrySpecification:
    from orpheus.derivations.common.xs_library import get_mixture

    geometry = StructuredGeometry.from_thicknesses(
        coord=CoordSystem.SPHERICAL, thicknesses=(0.5, 1.0, 0.5), mat_ids=(0, 1, 0), boundaries=(BC.reflective,))
    return GeometrySpecification(Materials({0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}), geometry, Eigen(_K))


_REFUSALS = [
    ("infinite-medium", lambda: InfiniteMediumSpecification(0, isotropic_mixture("A"), Eigen(_K)), "infinite"),
    ("fixed-source", lambda: GeometrySpecification(
        aba_specification(CoordSystem.SPHERICAL).materials, aba_specification(CoordSystem.SPHERICAL).geometry,
        FixedSource(RegionwiseConstant(np.ones((3, 2))))), "k"),
    ("eigen-along-scattering", lambda: GeometrySpecification(
        aba_specification(CoordSystem.SPHERICAL).materials, aba_specification(CoordSystem.SPHERICAL).geometry,
        Eigen(CellCoefficient.every(Channel.SCATTERING_EMISSION))), "k"),
    ("layered-slab", _slab_layered, "slab"),
    ("anisotropic-scattering", _anisotropic, "anisotropic"),
    ("homogeneous-sphere", lambda: GeometrySpecification(
        Materials({0: isotropic_mixture("A")}),
        StructuredGeometry(coord=CoordSystem.SPHERICAL, breakpoints=(0.0, 2.0), mat_ids=(0,), boundaries=(BC.reflective,)),
        Eigen(_K)), "layered"),
]


@pytest.mark.parametrize("build, fragment", [r[1:] for r in _REFUSALS], ids=[r[0] for r in _REFUSALS])
@pytest.mark.usefixtures("in_process")
def test_r7b2_2_refusals(build, fragment, monkeypatch: pytest.MonkeyPatch) -> None:
    """Refused at construction, before any solve. The anisotropic row: the trajectory resolvent solves isotropic
    scattering only, and ``get_mixture("B", "2g")`` carries a P1 moment (mean cosine 0.6); accepting it would
    silently solve another problem."""
    spy = _SolveSpy(monkeypatch)
    spec = build()
    with pytest.raises((TypeError, ValueError, NotImplementedError), match=f"(?i){fragment}"):
        api.reference(spec, coord=CoordSystem.SPHERICAL)
    assert spy.calls == 0


# ─────────────────────────────────────────────────────────────────────
# R7b2.3 — every reading uncertified; the verbs refuse before any solve
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("coord", _COORDS_SOLVING, ids=_IDS)
def test_r7b2_3_every_reading_is_uncertified(coord) -> None:
    ref = _aba(coord)
    observables = (
        Eigenvalue(),
        _group_total(0),
        FluxIntegral(RegionwiseConstant(np.array([[1.0, 0.0], [0.0, 0.0], [1.0, 0.0]]))),
        PointValue(0.8, 1),
        Ratio(_group_total(0), _group_total(1)),
    )
    for observable in observables:
        reading = ref.read(observable)
        assert type(reading) is Uncertified, (observable, reading)
        assert math.isfinite(reading.value) and reading.value != 0.0
    assert 1.0 < ref.read(Eigenvalue()).value < 1.6  # activation: an A|B|A k, not 1/k


class _Answer:
    def read(self, observable):
        return Measured(1.25)


@pytest.mark.parametrize("coord", _COORDS_SOLVING, ids=_IDS)
@pytest.mark.usefixtures("in_process")
def test_r7b2_3_the_verbs_refuse_before_any_solve(coord, monkeypatch: pytest.MonkeyPatch) -> None:
    spy = _SolveSpy(monkeypatch)
    ref = api.reference(aba_specification(coord))
    with pytest.raises(ReferenceNotValid, match="no certificate"):
        verify_agreement(_Answer(), Eigenvalue(), ref, 1e-3, NotYet(564, "no estimator"))
    with pytest.raises(ReferenceNotValid, match="no certificate"):
        verify_order([(1.0, _Answer()), (0.5, _Answer())], Eigenvalue(), ref, 1e-3, [NotYet(564, "x")] * 2, 2.0, 0.1)
    assert spy.calls == 0
    comparison = compare_uncertified(_Answer(), Eigenvalue(), ref, 1.0)
    assert comparison.reference_reading == ref.read(Eigenvalue())


# ─────────────────────────────────────────────────────────────────────
# R7b2.4 — the eigenvalue threads to the solver
# ─────────────────────────────────────────────────────────────────────


def _direct_k(coord: CoordSystem) -> float:
    """The multi-region solver called directly with the A|B|A arrays and the gates' parameters."""
    from tests.gates.sn.verification.analytical._aba_reference import ABA_MATERIAL_IDS, ABA_RADII, aba_materials

    mats = [aba_materials()[m] for m in ABA_MATERIAL_IDS]
    sigma_t = np.stack([np.asarray(m.SigT) for m in mats])
    sigma_s = np.stack([m.SigS[0].toarray() for m in mats])
    nu_sigma_f = np.stack([np.asarray(m.SigP) for m in mats])
    chi = np.stack([np.asarray(m.chi) for m in mats])
    module_name, name = api.SOLVERS[0 if coord is CoordSystem.SPHERICAL else 1]
    solver = getattr(importlib.import_module(module_name), name)
    result = solver(radii=np.array(ABA_RADII), sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f, chi=chi,
                    alpha=1.0, **api.QUADRATURE[coord], max_iter=api.MAX_ITER, tol=api.TOL, initial_k=api.INITIAL_K)
    return float(result.k_eff)


@pytest.mark.parametrize("coord", _COORDS_SOLVING, ids=_IDS)
def test_r7b2_4_the_eigenvalue_is_the_solvers_bit_for_bit(coord) -> None:
    """The specification's materials, geometry and laws reach the solver as the direct call spells them, and the
    eigenvalue is read off the result unaltered (the k chart)."""
    assert _aba(coord).read(Eigenvalue()).value == _direct_k(coord)


# ─────────────────────────────────────────────────────────────────────
# R7b2.5 — one transport
# ─────────────────────────────────────────────────────────────────────


def _oracle_and_profile(coord: CoordSystem):
    """A multi-region oracle on the A|B|A nodes and a non-trivial emission profile on them."""
    from orpheus.derivations.continuous.trajectory_resolvent import chord_oracle as co
    from orpheus.derivations.continuous.trajectory_resolvent.greens_function import _composite_per_region_gl
    from tests.gates.sn.verification.analytical._aba_reference import ABA_RADII

    radii = np.array(ABA_RADII)
    r_nodes, _, region = _composite_per_region_gl(radii, 12)
    sigma_t = np.array([0.5, 1.2, 0.5])
    profile = 1.0 + np.sin(3.0 * r_nodes) + 0.5 * region
    if coord is CoordSystem.SPHERICAL:
        mu = np.polynomial.legendre.leggauss(8)[0]
        oracle = co.MultiRegionSphereChordOracle(r_nodes, mu, 2.0, radii, sigma_t, 1.0, region)
    else:
        mu = np.polynomial.legendre.leggauss(4)[0]
        az = np.pi * (np.polynomial.legendre.leggauss(8)[0] + 1.0)
        oracle = co.MultiRegionCylinderChordOracle(r_nodes, mu, az, 2.0, radii, sigma_t, 1.0, region)
    return oracle, profile, r_nodes


@pytest.mark.rests_on(f"{_FIRST_LEG}[sphere]")
@pytest.mark.parametrize("coord", _COORDS, ids=_IDS)
def test_r7b2_5_the_oracle_at_its_nodes_is_its_default(coord) -> None:
    """``at=None`` (the solver's call) and ``at=r_nodes`` are one evaluation, bit for bit."""
    oracle, profile, r_nodes = _oracle_and_profile(coord)
    np.testing.assert_array_equal(api.apply_at(oracle, profile, None), api.apply_at(oracle, profile, r_nodes))


@pytest.mark.parametrize("coord", _COORDS, ids=_IDS)
def test_r7b2_5_the_knots_do_not_move_with_the_evaluation_points(coord) -> None:
    """Evaluating at a subset of the nodes gives exactly those rows of the full evaluation: the spline is built on
    the knots, never on the evaluation radii (the naive reuse of ``r_nodes`` for both reddens here)."""
    oracle, profile, r_nodes = _oracle_and_profile(coord)
    full = api.apply_at(oracle, profile, None)
    subset = np.arange(0, len(r_nodes), 3)
    np.testing.assert_array_equal(api.apply_at(oracle, profile, r_nodes[subset]), full[subset])


@pytest.mark.parametrize("coord", _COORDS_SOLVING, ids=_IDS)
def test_r7b2_5_the_uniform_body_reads_flat_at_the_infinite_medium_ratio(coord) -> None:
    """THEOREM (V_α1 on a layered uniform body): the reading is flat in r, off every node (0.13, 0.77, 1.31,
    1.99 cm), and the group ratio is the infinite medium's, exact. A transposed scattering matrix, a dropped 1/k
    or a group mix-up in the emission density moves the ratio (k∞ of fuel A is not 1)."""
    ref = _uniform(coord)
    exact = exact_infinite_medium_of(isotropic_mixture("A"))
    exact_ratio = float(Fraction(exact.flux[0]) / Fraction(exact.flux[1]))
    radii = (0.13, 0.77, 1.31, 1.99)
    for g in range(2):
        values = [ref.read(PointValue(r, g)).value for r in radii]
        assert max(values) - min(values) <= 1e-9 * abs(values[0]), (g, values)
    ratio = ref.read(Ratio(PointValue(0.77, 0), PointValue(0.77, 1))).value
    assert abs(ratio / exact_ratio - 1.0) <= 1e-9, (ratio, exact_ratio)
    assert abs(ref.read(Eigenvalue()).value / float(exact.k_inf) - 1.0) <= 1e-9


# ─────────────────────────────────────────────────────────────────────
# R7b2.6 — the reading's quadrature splits at the weight's steps
# ─────────────────────────────────────────────────────────────────────


_MEASURE = {
    CoordSystem.SPHERICAL: lambda a, b: 4.0 / 3.0 * math.pi * (b**3 - a**3),
    CoordSystem.CYLINDRICAL: lambda a, b: math.pi * (b * b - a * a),
}


def test_r7b2_6_symbolic_steps() -> None:
    """``Symbolic.steps``: the step locations inside the range, sorted; none for a smooth weight; a step on a
    non-polynomial argument refused with production's scope-edge message (one definition)."""
    r = Symbolic.r
    box = Symbolic.of(sympy.Piecewise((1, (r >= sympy.Rational(3, 5)) & (r < sympy.Rational(7, 10))), (0, True)))
    assert api.steps(box, 0.0, 2.0) == pytest.approx((0.6, 0.7), abs=0, rel=1e-15)
    assert api.steps(box, 0.0, 0.65) == pytest.approx((0.6,), abs=0, rel=1e-15)
    assert api.steps(Symbolic.of(sympy.Piecewise((1, r**2 < 2), (0, True))), 0.0, 2.0) == pytest.approx((math.sqrt(2.0),), rel=1e-15)
    assert api.steps(Symbolic.of(r**2 + 1), 0.0, 2.0) == ()
    with pytest.raises(ValueError, match="polynomial in r"):
        api.steps(Symbolic.of(sympy.Piecewise((1, sympy.sin(3 * r) > 0), (0, True))), 0.0, 2.0)


@pytest.mark.rests_on(f"{_HERE}::test_r7b2_6_symbolic_steps", f"{_HERE}::test_r7b2_5_the_uniform_body_reads_flat_at_the_infinite_medium_ratio")
@pytest.mark.parametrize("coord", _COORDS_SOLVING, ids=_IDS)
def test_r7b2_6_an_indicator_reads_its_measure_fraction(coord) -> None:
    """On the flat body, the indicator of (0.6, 0.7) cm (inside region 0, off every node) over the whole body
    reads the measure fraction exactly, to 1e-9; a quadrature that ignored the steps errs at first order in its
    panel width."""
    ref = _uniform(coord)
    measure = _MEASURE[coord]
    reading = ref.read(Ratio(_indicator(0.6, 0.7, 0), _group_total(0))).value
    assert abs(reading / (measure(0.6, 0.7) / measure(0.0, 2.0)) - 1.0) <= 1e-9, reading


@pytest.mark.parametrize("coord", _COORDS_SOLVING, ids=_IDS)
def test_r7b2_6_the_reading_is_additive_over_a_split(coord) -> None:
    """On the heterogeneous body, the indicator of (0.6, 0.7) reads the sum of (0.6, 0.65) and (0.65, 0.7)."""
    ref = _aba(coord)
    whole = ref.read(_indicator(0.6, 0.7, 1)).value
    parts = ref.read(_indicator(0.6, 0.65, 1)).value + ref.read(_indicator(0.65, 0.7, 1)).value
    assert abs(parts / whole - 1.0) <= 1e-9, (whole, parts)


@pytest.mark.usefixtures("in_process")
def test_r7b2_6_the_reading_consults_symbolic_steps(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE: the derivation finds a weight's steps through ``Symbolic.steps``, the one definition."""
    ref = _aba(CoordSystem.SPHERICAL)
    ref.read(Eigenvalue())  # solve first, outside the spy
    calls: list[tuple] = []
    original = getattr(Symbolic, "steps")  # the declared API (absent until the build lands)

    def counting(self, *args, **kwargs):
        calls.append(args)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Symbolic, "steps", counting)
    ref.read(_indicator(0.6, 0.7, 0))
    assert calls, "the reading never consulted Symbolic.steps"


# ─────────────────────────────────────────────────────────────────────
# R7b2.7 — the layers
# ─────────────────────────────────────────────────────────────────────


_FORBIDDEN = ("orpheus.transport", "orpheus.sn", "orpheus.diffusion", "orpheus.homogeneous", "orpheus.cp",
              "orpheus.moc", "orpheus.mc", "orpheus.kinetics", "orpheus.fuel", "orpheus.thermal_hydraulics",
              "orpheus.plotting")


def _imports(module_name: str) -> set[str]:
    path = pathlib.Path(importlib.util.find_spec(module_name).origin)  # type: ignore[union-attr]
    names: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            names |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.add(node.module)
    return names


def test_r7b2_7_the_reference_module_imports_no_higher_layer() -> None:
    """The new module (L0) reads the reference package and nothing at L2 or L3; positive control: it does import
    ``orpheus.reference``, so the AST walk saw its imports."""
    names = _imports(api.MODULE)
    assert any(n.startswith("orpheus.reference") for n in names), names
    assert not [n for n in names if n.startswith(_FORBIDDEN)], names


def test_r7b2_7_the_reference_package_imports_no_derivation() -> None:
    """``orpheus.reference`` stays below ``derivations`` (R6.7's law, re-read over every module of the package)."""
    import pkgutil

    import orpheus.reference as package

    seen = 0
    for info in pkgutil.iter_modules(package.__path__):
        names = _imports(f"orpheus.reference.{info.name}")
        seen += 1
        assert not [n for n in names if n.startswith("orpheus.derivations")], (info.name, names)
    assert seen >= 5


# ─────────────────────────────────────────────────────────────────────
# R7b2.2.1–R7b2.2.4 — the review round of 7b.2.2 (qa F1–F3, elegance B1)
# ─────────────────────────────────────────────────────────────────────
#
# Ids R7b2.2.<n>: R7b2.8 and R7b2.9 already name the SN reading's gates, so the review round extends R7b2.2's family.


_OFF_NODE = (0.31, 1.2, 1.73)


def _brute_sphere(ref, group: int, r: float, n_mu: int = 2000) -> float:
    """φ_g(r) through the reference's OWN oracle and emission density, with an UNSPLIT Gauss–Legendre μ rule of
    ``n_mu`` points on [-1, 1] (2π ∫ ψ dμ): independent of the reading's split angular rule and of its fibre."""
    import dataclasses

    mu, w = np.polynomial.legendre.leggauss(n_mu)
    rays = dataclasses.replace(api.oracle(ref, group), mu_nodes=mu)
    psi = rays.apply_operator(api.emission_density(ref)[group], 0.0, n_traj_quad=64, at=np.array([r]))
    return float(2.0 * np.pi * (psi[0] @ w))


def _brute_cylinder(ref, group: int, r: float, n_theta: int, n_azimuth: int) -> float:
    """φ_g(r) = ∫ sinθ dθ ∫ ψ dφ by UNSPLIT Gauss–Legendre in θ ∈ [0, π] and φ ∈ [0, 2π], the reference's oracle."""
    import dataclasses

    x, wx = np.polynomial.legendre.leggauss(n_theta)
    theta, w_theta = 0.5 * np.pi * (x + 1.0), 0.5 * np.pi * wx
    y, wy = np.polynomial.legendre.leggauss(n_azimuth)
    azimuth, w_azimuth = np.pi * (y + 1.0), np.pi * wy
    rays = dataclasses.replace(api.oracle(ref, group), mu_axial_nodes=np.cos(theta), phi_az_nodes=azimuth)
    psi = rays.apply_operator(api.emission_density(ref)[group], 0.0, n_traj_quad=64, at=np.array([r]))
    return float((w_theta * np.sin(theta)) @ psi[0] @ w_azimuth)


def test_r7b2_2_1_the_sphere_reading_against_an_unsplit_fine_angular_rule() -> None:
    """qa F1 (a), VALUE on the heterogeneous body: the reading at three off-node radii, both groups, against an
    unsplit 2000-point μ rule through the same oracle, to 1e-5 relative (``[M]`` 2026-10-03: the unsplit rule's own
    error dominates, ~n^-1.5 at r = 1.73: 9.6e-6 at 1000 points, 1.5e-6 at 3000, about 2.7e-6 at 2000). An
    angular measure off by a direction-dependent factor (integrating μ ∈ [0, 1] and doubling) moves these by up
    to 9 % on A|B|A; the uniform body's rows are blind to it (an isotropic ψ)."""
    ref = _aba(CoordSystem.SPHERICAL)
    for g in range(2):
        for r in _OFF_NODE:
            reading = ref.read(PointValue(r, g)).value
            brute = _brute_sphere(ref, g, r)
            assert abs(reading / brute - 1.0) <= 1e-5, (g, r, reading, brute)


@pytest.mark.slow
def test_r7b2_2_1_the_cylinder_reading_against_an_unsplit_fine_angular_rule() -> None:
    """qa F1 (a) on the cylinder, at the gates' resolution: dropping sinθ moves the A|B|A ratios by up to 8.7 %."""
    ref = _aba(CoordSystem.CYLINDRICAL)
    for g in range(2):
        for r in _OFF_NODE:
            reading = ref.read(PointValue(r, g)).value
            brute = _brute_cylinder(ref, g, r, n_theta=96, n_azimuth=768)
            assert abs(reading / brute - 1.0) <= 1e-5, (g, r, reading, brute)


@functools.cache
def _cheap_cylinder():
    """A|B|A cylinder at (4, 2, 4) with a cheap reading (4 points per angular piece, 16 per ray segment)."""
    ref = api.reference(aba_specification(CoordSystem.CYLINDRICAL),
                        quadrature={"n_r": 4, "n_mu_axial": 2, "n_phi_az": 4, "n_traj_quad": 16})
    return api.with_reading(ref, angular_points_per_piece=4, ray_points_per_segment=16)


def test_r7b2_2_1_a_cheap_cylinder_reading_law_runs_in_not_slow() -> None:
    """qa F2: one cylinder READING law in ``not slow`` (``_CylinderRays.scalar_flux`` ran 0 times there): at one
    off-node radius, group 1, the cheap reading against an unsplit (32, 256) rule. ``[M]`` 2026-10-03: they differ by
    7.4e-5, the cheap reading's own angular error (the brute moves 1e-6 to (96, 768)); the band 1e-3 is about ten
    times that, and dropping sinθ (8.7 %) is far outside it."""
    ref = _cheap_cylinder()
    reading = ref.read(PointValue(1.2, 1)).value
    brute = _brute_cylinder(ref, 1, 1.2, n_theta=32, n_azimuth=256)
    assert abs(reading / brute - 1.0) <= 1e-3, (reading, brute)


@pytest.mark.parametrize("coord", _COORDS_SOLVING, ids=_IDS)
def test_r7b2_2_2_the_emission_density_is_the_solves_fixed_point(coord) -> None:
    """qa F1 (b), the fixed-point identity on the HETEROGENEOUS body: the oracle with the solve's own knots and
    angular rule, applied to the reference's emission density (the solve's LAST source, hoisted onto its result),
    reproduces the solve's final ψ up to one scalar (the gauge, the declared stabiliser), with a residual at the
    rounding level, 1e-12 relative to max |ψ|. It pins the density (σ_s orientation, χ, 1/k, the group order)
    where the uniform body cannot, and that the density read is the one the solve used. Structurally blind to the
    READING's angular rule, which it never calls (R7b2.2.1 is that rule's gate)."""
    ref = _aba(coord)
    module_name, name = api.SOLVERS[0 if coord is CoordSystem.SPHERICAL else 1]
    from tests.gates.sn.verification.analytical._aba_reference import ABA_MATERIAL_IDS, ABA_RADII, aba_materials

    mats = [aba_materials()[m] for m in ABA_MATERIAL_IDS]
    raw = getattr(importlib.import_module(module_name), name)(
        radii=np.array(ABA_RADII), sigma_t=np.stack([np.asarray(m.SigT) for m in mats]),
        sigma_s=np.stack([m.SigS[0].toarray() for m in mats]), nu_sigma_f=np.stack([np.asarray(m.SigP) for m in mats]),
        chi=np.stack([np.asarray(m.chi) for m in mats]), alpha=1.0, **api.QUADRATURE[coord],
        max_iter=api.MAX_ITER, tol=api.TOL, initial_k=api.INITIAL_K,
    )
    density = api.emission_density(ref)
    transported = np.stack([
        api.apply_at(api.oracle(ref, g), density[g], None, n_traj_quad=api.QUADRATURE[coord]["n_traj_quad"])
        for g in range(2)
    ])
    psi = np.asarray(raw.psi_g)
    scale = float(np.vdot(transported, psi) / np.vdot(psi, psi))
    residual = float(np.max(np.abs(transported - scale * psi))) / float(np.max(np.abs(transported)))
    assert residual <= 1e-12, (residual, scale)


def test_r7b2_2_3_steps_at_algebraic_and_transcendental_locations() -> None:
    """qa F3: a step at π/4 and at √2/2 (SymPy's polynomial root finder over QQ[π] or EX raised on these)."""
    r = Symbolic.r
    quarter_pi = Symbolic.of(sympy.Piecewise((1, r < sympy.pi / 4), (0, True)))
    root_half = Symbolic.of(sympy.Piecewise((1, r < sympy.sqrt(2) / 2), (0, True)))
    assert api.steps(quarter_pi, 0.0, 2.0) == pytest.approx((math.pi / 4,), rel=1e-15)
    assert api.steps(root_half, 0.0, 2.0) == pytest.approx((math.sqrt(2.0) / 2,), rel=1e-15)


def test_r7b2_2_3_the_mesh_integrates_a_step_at_pi_over_4() -> None:
    """qa F3, the regression on ``Mesh1D.cell_integrals``: ``r < π/4`` on the A|B|A sphere's mesh totals
    4π/3·(π/4)³ (2.029356063208384 at ``fa38de31``), to 64 ulp."""
    from orpheus.mesh import CellsByCount, Mesher
    from tests.gates.sn.verification.analytical._aba_reference import aba_geometry

    mesh = Mesher(aba_geometry(CoordSystem.SPHERICAL)).partition(tuple(CellsByCount.uniform_width(n) for n in (2, 4, 2))).mesh
    weight = Symbolic.of(sympy.Piecewise((1, Symbolic.r < sympy.pi / 4), (0, True)))
    total = float(np.sum(mesh.cell_integrals(weight)))
    expected = 4.0 / 3.0 * math.pi * (math.pi / 4) ** 3
    assert abs(total - expected) <= 64 * math.ulp(expected), (total, expected)


@pytest.mark.parametrize("make", [
    lambda r: sympy.arg(r - 1),
    lambda r: sympy.atan2(r - 1, 0),
], ids=["arg", "atan2"])
def test_r7b2_2_3_an_unlisted_non_smooth_construct_is_refused(make) -> None:
    """``Symbolic.steps`` is an ALLOW-LIST: a construct whose jumps it does not locate (``arg(r − 1)`` jumps at 1,
    ``atan2(r − 1, 0)`` too) is refused, never integrated across silently."""
    with pytest.raises(ValueError, match=api.UNLOCATED_REFUSAL):
        api.steps(Symbolic.of(make(Symbolic.r)), 0.0, 2.0)


def test_r7b2_2_4_the_public_constructor_derives_its_rays_from_the_billiard() -> None:
    """Elegance B1: the derivation's public constructor takes the specification and derives its billiard and its
    rays from it (no rays passed beside them, so a sphere billiard cannot be read with the cylinder's rays). Built
    directly, it reads what the factory's reference reads, bit for bit."""
    from orpheus.reference.solution import ReferenceSolution

    spec = aba_specification(CoordSystem.SPHERICAL)
    built = api.derivation(spec)
    init_fields = {f.name for f in __import__("dataclasses").fields(type(built)) if f.init}
    assert "rays" not in init_fields and "billiard" not in init_fields, init_fields
    direct = ReferenceSolution(spec, built, None)
    assert direct.read(Eigenvalue()).value == _aba(CoordSystem.SPHERICAL).read(Eigenvalue()).value
    assert direct.read(PointValue(1.2, 0)).value == _aba(CoordSystem.SPHERICAL).read(PointValue(1.2, 0)).value


@pytest.mark.parametrize("build, fragment", [r[1:] for r in _REFUSALS], ids=[r[0] for r in _REFUSALS])
@pytest.mark.usefixtures("in_process")
def test_r7b2_2_4_the_public_constructor_refuses_what_the_factory_refuses(build, fragment, monkeypatch) -> None:
    """The constructor is the one door: every refusal of R7b2.2 holds when the derivation is built directly."""
    spy = _SolveSpy(monkeypatch)
    with pytest.raises((TypeError, ValueError, NotImplementedError), match=f"(?i){fragment}"):
        api.derivation(build(), CoordSystem.SPHERICAL)
    assert spy.calls == 0


# ─────────────────────────────────────────────────────────────────────
# Content identity: the derivation is the key of its memoised readings (#405 P3)
# ─────────────────────────────────────────────────────────────────────


_ROSTER_QUADRATURE = {"n_r": 8, "n_mu": 8, "n_traj_quad": 16}


def _aba_sphere(thicknesses: tuple[float, ...] = (0.5, 1.0, 0.5)) -> GeometrySpecification:
    """The A|B|A sphere's specification, built FRESH on every call (S5.6 rebuilds each subject under a decoy
    encoder; the cached ``aba_specification`` would carry digests from the honest one)."""
    from tests.gates.sn.verification.analytical._aba_reference import ABA_MATERIAL_IDS

    materials = Materials({0: isotropic_mixture("A"), 1: isotropic_mixture("B")})
    geometry = StructuredGeometry.from_thicknesses(
        coord=CoordSystem.SPHERICAL, thicknesses=thicknesses, mat_ids=ABA_MATERIAL_IDS, boundaries=(BC.reflective,),
    )
    return GeometrySpecification(materials, geometry, Eigen(CellCoefficient.every(Channel.FISSION_EMISSION)))


def _roster_derivation(**changes: Any) -> Any:
    """The sphere A|B|A derivation at the gates' resolution, with ``changes`` to its init fields."""
    fields = {"specification": _aba_sphere(), "solver_quadrature": _ROSTER_QUADRATURE,
              "max_iter": 500, "tol": 1e-10, "initial_k": 1.0} | changes
    return api.derivation_class()(**fields)


def _reading_quadrature(**changes: int) -> Any:
    return api.module().ReadingQuadrature(**changes)


ROSTER: tuple[Entry, ...] = (
    Entry(
        cls=api.derivation_class(),
        base=_roster_derivation,
        parts=("specification", "solver_quadrature", "max_iter", "tol", "initial_k", "quadrature"),
        perturb={
            "specification": (leg("another body", lambda: _roster_derivation(specification=_aba_sphere((0.5, 1.0, 0.55)))),),
            "solver_quadrature": (leg("finer radial nodes", lambda: _roster_derivation(solver_quadrature={**_ROSTER_QUADRATURE, "n_r": 10})),),
            "max_iter": (leg("another budget", lambda: _roster_derivation(max_iter=400)),),
            "tol": (leg("another tolerance", lambda: _roster_derivation(tol=1e-9)),),
            "initial_k": (leg("another guess", lambda: _roster_derivation(initial_k=1.2)),),
            "quadrature": (leg("more ray points", lambda: _roster_derivation(quadrature=_reading_quadrature(ray_points_per_segment=32))),),
        },
        pairs=(("a-mapping-and-its-frozen-twin", _roster_derivation, lambda: _roster_derivation(solver_quadrature=dict(_ROSTER_QUADRATURE))),),
    ),
)


@pytest.mark.parametrize("entry", ROSTER, ids=[e.id for e in ROSTER])
def test_p3_the_derivation_population_is_its_parts(entry: Entry) -> None:
    """The derivation's content is its six init fields; the ``Billiard`` and the rays class are derived."""
    check_population(entry)


@pytest.mark.parametrize("entry, part, the_leg", perturbation_ids(ROSTER),
                         ids=[param_id(e.id, p, lg[0]) for e, p, lg in perturbation_ids(ROSTER)])
def test_p3_each_init_field_moves_the_derivation_digest(entry: Entry, part: str, the_leg: Any) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.parametrize("entry, pair", pair_ids(ROSTER), ids=[param_id(e.id, p[0]) for e, p in pair_ids(ROSTER)])
def test_p3_equal_fields_are_one_derivation(entry: Entry, pair: Any) -> None:
    check_equal_pair(entry, pair)
