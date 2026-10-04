r"""The SN production reading: ``Solution.read(observable) -> Measured`` (#405 P2 step 7b.2.1, gates R7b2.9).

The forward S\ :sub:`N` solution reads an observable as production's
self-report, the twin of ``HomogeneousResult.read`` (G1: the answer is the
receiver of ``read``):

* ``Eigenvalue`` reads ``outcome.keff`` on an ``EigenOutcome``; a
  ``SourceOutcome`` refuses it;
* ``FluxIntegral(w)`` reads :math:`\sum_{g,i}\phi_{g,i}\,
  \texttt{mesh.cell\_integrals}(w)_{g,i}`, φ the cell averages
  ``scalar_flux.values`` in the solution's own gauge (exact for a
  piecewise-constant answer);
* ``Ratio`` reads ``ratio.quotient(self.read)`` (the one ratio rule, R7.19);
* ``PointValue`` is refused (a cell-average answer has no point value), and
  a 2-D problem refuses a flux integral (#569); ``AdjointSolution`` has no
  ``read`` in this step.

Gate ids ``R7b2.9.<n>`` (``.claude/plans/reference_p2_spec.md`` §1.7b.2,
"7b.2.1"). First red on ``fa38de31``: ``AttributeError``, ``Solution`` has no
``read``. Fixture: the A|B|A sphere, 2 groups, (2, 4, 2) cells, Gauss–Legendre
4, about 1 s. Claim kind THEOREM unless marked; every row ``foundation``
except the end-to-end reduction row (``l1``).
"""
from __future__ import annotations

import functools
import math
import warnings

import numpy as np
import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.derivations.common.exact_homogeneous import exact_infinite_medium_reference
from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic
from orpheus.numerics.observable import Eigenvalue, FluxIntegral, PointValue, Ratio
from orpheus.numerics.outcome import Measured
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.question import Eigen
from orpheus.sn import solve_sn, solve_sn_fixed_source
from orpheus.sn.solution import AdjointSolution, Solution
from orpheus.specification.specification import InfiniteMediumSpecification

_HERE = "tests/gates/sn/test_solution_read.py"
_CELLS = "tests/gates/mesh/test_mesh1d_cell_integrals.py"
_RATIO_RULE = "tests/gates/reference/test_verification.py::test_r7_19_every_answer_reads_a_ratio_through_the_one_quotient"
_REGIONS = (0.5, 1.5, 2.0)


@functools.cache
def _aba():
    """The A|B|A sphere eigen solve and its mesh."""
    g = StructuredGeometry.from_thicknesses(
        coord=CoordSystem.SPHERICAL, thicknesses=(0.5, 1.0, 0.5), mat_ids=(0, 1, 0), boundaries=(BC.reflective,),
    )
    mesh = Mesher(g).partition(tuple(CellsByCount.uniform_width(n) for n in (2, 4, 2))).mesh
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = solve_sn({0: get_mixture("A", "2g"), 1: get_mixture("B", "2g")}, mesh,
                       Quadrature.gauss_legendre(n_ordinates=4), keff_tol=1e-10, flux_tol=1e-9,
                       max_inner=200, inner_tol=1e-11)
    return sol, mesh


def _phi_volumes():
    sol, mesh = _aba()
    return np.asarray(sol.scalar_flux.values, dtype=float), np.asarray(mesh.volumes, dtype=float)


def _group_total(g: int) -> FluxIntegral:
    return FluxIntegral(Symbolic.of(*(1 if h == g else 0 for h in range(2))))


def _cell_indicator(i: int, g: int) -> FluxIntegral:
    """Group ``g``'s indicator of cell ``i`` (half-open), unnormalised."""
    _, mesh = _aba()
    a, b = (float(x) for x in np.asarray(mesh.edges)[i:i + 2])
    w = Symbolic.r
    import sympy
    step = sympy.Piecewise((1, (w >= a) & (w < b)), (0, True))
    return FluxIntegral(Symbolic.of(*(step if h == g else 0 for h in range(2))))


# ─────────────────────────────────────────────────────────────────────
# R7b2.9.1 — the eigenvalue
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_r7b2_9_1_k_is_the_outcome_k_bit_for_bit() -> None:
    sol, _ = _aba()
    reading = sol.read(Eigenvalue())
    assert type(reading) is Measured
    assert reading.value == sol.outcome.keff  # the k chart, not 1/k; no re-computation


@pytest.mark.foundation
def test_r7b2_9_1_a_source_answer_refuses_the_eigenvalue() -> None:
    slab = Mesher(StructuredGeometry.slab((0.0, 1.0), (0,), left=BC.reflective, right=BC.vacuum)).partition(
        CellsByCount.uniform_width(4)).mesh
    quadrature = Quadrature.gauss_legendre(n_ordinates=4)
    fs = solve_sn_fixed_source({0: get_mixture("B", "2g")}, slab, quadrature, np.ones((4, 2, 4)),
                               max_inner=400, inner_tol=1e-12)
    with pytest.raises(ValueError, match="eigen"):
        fs.read(Eigenvalue())
    assert fs.read(_group_total(0)).value > 0.0  # the source answer still reads a flux integral


# ─────────────────────────────────────────────────────────────────────
# R7b2.9.2 — flux integrals
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_CELLS}::test_r7b2_8_1_the_layout_and_the_indicator_of_each_cell")
def test_r7b2_9_2_a_group_total_is_the_volume_weighted_sum() -> None:
    """Group g's weight 1 reads Σ_i φ_{g,i} V_i (the cell integrals of 1 are the volumes to a few ulp)."""
    sol, _ = _aba()
    phi, volumes = _phi_volumes()
    for g in range(2):
        expected = float(np.sum(phi[g] * volumes))
        got = sol.read(_group_total(g)).value
        assert abs(got - expected) <= 16 * len(volumes) * math.ulp(expected), (g, got, expected)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_HERE}::test_r7b2_9_2_a_group_total_is_the_volume_weighted_sum")
def test_r7b2_9_2_the_cells_sum_to_the_whole_domain() -> None:
    """Σ over every cell and group of the per-cell indicator readings equals the whole-domain reading."""
    sol, mesh = _aba()
    parts = sum(sol.read(_cell_indicator(i, g)).value for i in range(mesh.N) for g in range(2))
    whole = sol.read(FluxIntegral(Symbolic.of(1, 1))).value
    assert abs(parts - whole) <= 64 * mesh.N * math.ulp(whole), (parts, whole)


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_CELLS}::test_r7b2_8_5_a_region_table_reads_through_the_labels")
def test_r7b2_9_2_a_region_table_reads_regions_0_and_2_separately() -> None:
    """The indicator of region 0 and the indicator of region 2 (both material 0) read the flux over each region
    alone; the reference finds each cell's region by containment of its centre."""
    sol, mesh = _aba()
    phi, volumes = _phi_volumes()
    centres = 0.5 * (np.asarray(mesh.edges)[:-1] + np.asarray(mesh.edges)[1:])
    region = np.searchsorted(_REGIONS, centres)
    for k in (0, 2):
        table = np.zeros((3, 2))
        table[k, 0] = 1.0
        expected = float(np.sum((phi[0] * volumes)[region == k]))
        got = sol.read(FluxIntegral(RegionwiseConstant(table))).value
        assert abs(got - expected) <= 16 * mesh.N * math.ulp(expected), (k, got, expected)
    r0 = sol.read(FluxIntegral(RegionwiseConstant(np.array([[1.0, 0.0], [0.0, 0.0], [0.0, 0.0]])))).value
    r2 = sol.read(FluxIntegral(RegionwiseConstant(np.array([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0]])))).value
    assert r0 != r2, "regions 0 and 2 read as one"


@pytest.mark.foundation
@pytest.mark.rests_on(f"{_HERE}::test_r7b2_9_2_a_group_total_is_the_volume_weighted_sum")
def test_r7b2_9_2_the_reading_reads_the_flux(monkeypatch: pytest.MonkeyPatch) -> None:
    """X1 negative: φ scaled by 1 + 1e-11 moves the flux-integral reading by about 1e-11 relative, and leaves k."""
    sol, _ = _aba()
    obs = _group_total(1)
    honest, k = sol.read(obs).value, sol.read(Eigenvalue()).value
    original = type(sol).scalar_flux

    class _Scaled:
        def __init__(self, inner):
            self.values = np.asarray(inner.values) * (1.0 + 1e-11)

    monkeypatch.setattr(type(sol), "scalar_flux", property(lambda self: _Scaled(original.__get__(self))))
    moved = sol.read(obs).value
    assert abs(moved / honest - 1.0 - 1e-11) <= 1e-13, (honest, moved)
    assert sol.read(Eigenvalue()).value == k


# ─────────────────────────────────────────────────────────────────────
# R7b2.9.3 — ratios, through the one rule; group order
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_RATIO_RULE)
def test_r7b2_9_3_a_ratio_reads_through_ratio_quotient(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE: rebinding ``Ratio.quotient`` to divide the other way moves the SN reading to its reciprocal, the call
    coming from the SN solution: it keeps no private ratio rule (R7.19 extended to this answer)."""
    sol, _ = _aba()
    ratio = Ratio(_group_total(0), _group_total(1))
    honest = sol.read(ratio).value
    calls: list[str] = []

    def reversed_quotient(self, read):
        calls.append(type(read.__self__).__name__)
        return read(self.denominator) / read(self.numerator)

    monkeypatch.setattr(Ratio, "quotient", reversed_quotient)
    moved = sol.read(ratio).value
    assert calls == ["Solution"], calls
    assert abs(moved * honest - 1.0) <= 1e-14 and honest != 1.0


@pytest.mark.foundation
def test_r7b2_9_3_the_group_order() -> None:
    """On the 2-group fixture the fast/thermal ratio is Σφ₀V / Σφ₁V, not its reciprocal (≠ 1, so a group
    reversal anywhere on the path reddens)."""
    sol, _ = _aba()
    phi, volumes = _phi_volumes()
    expected = float(np.sum(phi[0] * volumes)) / float(np.sum(phi[1] * volumes))
    got = sol.read(Ratio(_group_total(0), _group_total(1))).value
    assert abs(expected - 1.0) > 0.1
    assert abs(got / expected - 1.0) <= 1e-13, (got, expected)


# ─────────────────────────────────────────────────────────────────────
# R7b2.9.4 — refusals and scope
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_r7b2_9_4_a_point_value_is_refused() -> None:
    sol, _ = _aba()
    with pytest.raises(ValueError, match="point value"):
        sol.read(PointValue(1.0, 0))


@pytest.mark.foundation
def test_r7b2_9_4_a_weight_of_another_group_count_is_refused() -> None:
    """A weight's group count is the answer's: a 3-group weight fails loudly only by accident (a broadcast error),
    and a 1-group weight BROADCASTS SILENTLY over both groups (``[M]`` 2026-10-03 on the first code: it read
    12.87859622934138, the two-group total). Both are refused, keyed."""
    sol, _ = _aba()
    for weight in (Symbolic.of(1), Symbolic.of(1, 1, 1), RegionwiseConstant(np.ones((3, 1)))):
        with pytest.raises(ValueError, match="group"):
            sol.read(FluxIntegral(weight))


@pytest.mark.foundation
def test_r7b2_9_4_a_2d_problem_refuses_a_flux_integral() -> None:
    """The 2-D cell integrals and labels are #569: refused, the issue named."""
    mesh = Mesh2D(
        edges_x=np.array([0.0, 0.5, 1.0]), edges_y=np.array([0.0, 0.5, 1.0]), mat_map=np.zeros((2, 2), dtype=int),
        coord=CoordSystem.CARTESIAN,
        face_laws={"xmin": BC("reflective"), "xmax": BC("vacuum"), "ymin": BC("reflective"), "ymax": BC("vacuum")},
    )
    quadrature = Quadrature.level_symmetric(sn_order=4)
    n = len(np.asarray(quadrature.weights))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fs = solve_sn_fixed_source({0: get_mixture("B", "2g")}, mesh, quadrature, np.ones((n, 2, 2, 2)),
                                   max_inner=200, inner_tol=1e-10, inner_schedule="jacobi")
    with pytest.raises((ValueError, NotImplementedError), match="#569"):
        fs.read(_group_total(0))


@pytest.mark.foundation
def test_r7b2_9_4_the_adjoint_has_no_read_yet() -> None:
    """Declared scope: the adjoint solution's reading (an importance, not a flux) is not built in this step."""
    assert hasattr(Solution, "read") and not hasattr(AdjointSolution, "read")


# ─────────────────────────────────────────────────────────────────────
# R7b2.9.5 — end to end: the homogeneous reduction (k = k∞), a theorem, not a verification
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_r7b2_9_5_a_reflective_homogeneous_slab_reads_k_infinity() -> None:
    """THEOREM (the reduction identity): a homogeneous body with reflective faces has the flat flux as its exact
    solution for any spatial and angular discretisation, so its k is k∞. The SN reading is compared with the
    exact infinite medium's reading at 10 × the solve's ``keff_tol``.

    Not ``verify_agreement``: the SN answer answers a slab's ``GeometrySpecification``, the reference an
    ``InfiniteMediumSpecification``, and equating the two is this theorem, a cross-specification pairing
    nothing checks until P4's projection (spec §1.7b.2, refuted candidates). ``[M]`` 2026-10-03: |Δk| = 6.8e-14
    at keff_tol 1e-12. Blind to every spatial and angular operator (flat flux): it pins the eigenvalue reading's
    chart and threading, not transport."""
    keff_tol = 1e-12
    mesh = Mesher(StructuredGeometry.slab((0.0, 3.0), (0,), left=BC.reflective, right=BC.reflective)).partition(
        CellsByCount.uniform_width(6)).mesh
    mixture = get_mixture("A", "2g")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = solve_sn({0: mixture}, mesh, Quadrature.gauss_legendre(n_ordinates=4), keff_tol=keff_tol,
                       flux_tol=1e-10, max_inner=200, inner_tol=1e-12)
    reference = exact_infinite_medium_reference(
        InfiniteMediumSpecification(0, mixture, Eigen(CellCoefficient.every(Channel.FISSION_EMISSION))),
    )
    k_ref = reference.read(Eigenvalue()).value
    k_sn = sol.read(Eigenvalue()).value
    assert abs(k_sn - k_ref) <= 10 * keff_tol * k_ref, (k_sn, k_ref)
    assert abs(1.0 / k_sn - k_ref) > 0.1  # activation: the reciprocal chart would be red


# ─────────────────────────────────────────────────────────────────────
# R7b2.9.6 — the review round: one 2-D route, and a problem with no mesh
# ─────────────────────────────────────────────────────────────────────


def _two_d_source_answer():
    mesh = Mesh2D(
        edges_x=np.array([0.0, 0.5, 1.0]), edges_y=np.array([0.0, 0.5, 1.0]), mat_map=np.zeros((2, 2), dtype=int),
        coord=CoordSystem.CARTESIAN,
        face_laws={"xmin": BC("reflective"), "xmax": BC("vacuum"), "ymin": BC("reflective"), "ymax": BC("vacuum")},
    )
    quadrature = Quadrature.level_symmetric(sn_order=4)
    n = len(np.asarray(quadrature.weights))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return solve_sn_fixed_source({0: get_mixture("B", "2g")}, mesh, quadrature, np.ones((n, 2, 2, 2)),
                                     max_inner=200, inner_tol=1e-10, inner_schedule="jacobi")


@pytest.mark.foundation
def test_r7b2_9_6_the_2d_refusal_is_the_meshs(monkeypatch: pytest.MonkeyPatch) -> None:
    """ROUTE: ``Solution.read`` keeps no 2-D branch of its own; it asks the mesh. Rebinding
    ``Mesh2D.cell_integrals`` to a decoy that returns the co-vector of ones makes the 2-D read SUCCEED with
    Σ φ, so the refusal lives in the mesh alone."""
    fs = _two_d_source_answer()
    phi = np.asarray(fs.scalar_flux.values, dtype=float)
    monkeypatch.setattr(Mesh2D, "cell_integrals", lambda self, weight: np.ones_like(phi))
    assert fs.read(_group_total(0)).value == pytest.approx(float(np.sum(phi)), rel=1e-15, abs=0.0)


@pytest.mark.foundation
def test_r7b2_9_6_a_problem_with_no_mesh_refuses_a_flux_integral() -> None:
    """A 3-D problem posed from axes holds no mesh (``SNProblem.mesh is None``): its eigenvalue reads, its flux
    integral is refused."""
    from orpheus.mesh import AxisMesh

    axes = tuple(AxisMesh(edges=np.linspace(0.0, ext, n + 1), bc_low=BC.reflective, bc_high=BC.reflective)
                 for ext, n in ((0.5, 2), (0.75, 3), (0.5, 2)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = solve_sn({0: get_mixture("A", "2g")}, axes, Quadrature.level_symmetric(sn_order=4),
                       keff_tol=1e-8, inner_tol=1e-9)
    assert sol.problem.mesh is None
    assert sol.read(Eigenvalue()).value == sol.outcome.keff
    with pytest.raises(ValueError, match="holds no mesh"):
        sol.read(_group_total(0))
