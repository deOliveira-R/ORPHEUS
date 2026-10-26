r"""CP's first-collision probability against closed forms, and CP's vacuum eigenvalue against its own pencil (#578, P2 group D).

**D1.** ``P_cell[i, j]`` is the probability that a neutron born uniformly and
isotropically in region i has its first collision in region j, with no
boundary closure. Summed over j and volume-averaged over a homogeneous body,
it is the body's collision probability ``P_c = 1 - P_esc``, whose closed
forms are written in ``tests/gates/derivations/_characteristic_mp.py`` from
the chord-length distributions (sphere: Hebert's
``3/(8 tau^3) (2 tau^2 - 1 + (1 + 2 tau) e^{-2 tau})``; cylinder ``(1 -
P_ss)/(2 tau)`` with Bickley's ``P_ss``; slab ``(1 - 2 E_3)/(2 tau)``). They
share no line integral, no impact-parameter rule and no density with CP.

Until #578 the sphere's density over the impact parameter was spelled ``y``
with a literal 2 (the cylinder's two lines at ``+-y``) where the measure is
``2 pi y`` (``Chart.beam_density``): production and the ``flat_source_cp``
reference both read ``1/pi`` of the closed form (`[M]` 2026-10-10, ratio
3.1416 at tau = 0.1, 1, 5), and agreed with each other because they share the
spelling (X4). White-boundary renormalisation hid it from every k gate; it was
#509's cause (``test_white_boundary_infinite_medium.py``).

**The vacuum eigenvalue.** ``CPSolver.compute_keff`` divides production by
the net removal with no leakage term, a premise true only under white, so its
k under vacuum is not the eigenvalue of the CP pencil the iteration solves.
The dense dominant eigenvalue of ``(I - T S)^{-1} T F``, built from the same
``P_inf`` matrices, is the reference (W2 memo
``scratch/characteristic_architecture/p2/cp_pi_investigation/memo.md``,
probes C, E, H). Not fixed by the user's ruling (CP's pencil is a later
campaign): a strict xfail, whose strictness is asserted by introspection.

The reference P_cell is read through a closure that never returns
(``surface_area`` -> 1e300, so ``P_in -> 0`` and ``P_inf = P_cell``);
production's under a vacuum outer law, where ``P_inf = P_cell``.
"""
from __future__ import annotations

import contextlib
import io

import numpy as np
import pytest

from orpheus.cp.solver import CPMesh, CPParams, solve_cp
from orpheus.data.cell_xs import assemble_cell_xs
from orpheus.data.macro_xs.synthetic import get_mixture
from orpheus.derivations.continuous.flat_source_cp import geometry as fsg
from orpheus.geometry import CoordSystem
from orpheus.mesh import CellsByCount, Mesher
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates.cp._cp_geometry import cp_body
from tests.gates.derivations import _characteristic_mp as closed
from tests.gates.geometry._measure import CONSTANT as _MEASURE_CONSTANT

_HERE = "tests/gates/cp/test_cp_first_collision.py::"
_CHART = {CoordSystem.SPHERICAL: "sphere", CoordSystem.CYLINDRICAL: "cylinder", CoordSystem.CARTESIAN: "slab"}
_KIND = {CoordSystem.SPHERICAL: "sphere-1d", CoordSystem.CYLINDRICAL: "cylinder-1d", CoordSystem.CARTESIAN: "slab"}


class _NoClosure(fsg.FlatSourceCPGeometry):
    """The reference geometry with a white closure that never returns: ``P_in = Sigma V P_out / S -> 0``."""

    def surface_area(self, R: float) -> float:
        return 1e300


def _truth(coord: CoordSystem, tau: float) -> float:
    """The closed-form ``P_c`` of the unit body at ``Sigma = tau``. CP's slab is half of a symmetric slab (a mirror
    at its left face, ``cp_body``), so its body is the slab of width 2: optical width ``2 tau``."""
    width = 2.0 * tau if coord is CoordSystem.CARTESIAN else tau
    return 1.0 - float(closed.escape_probability(_CHART[coord], width))


def _production_pc(coord: CoordSystem, tau: float) -> float:
    geom = cp_body(coord, (0.0, 0.5, 1.0), (0, 0), 0)
    mesh = Mesher(geom).partition((CellsByCount.uniform_volume(3), CellsByCount.uniform_volume(4))).mesh
    P = CPMesh(mesh, boundary_conditions={0: WallCondition.vacuum()}).compute_pinf_group(np.full(mesh.N, tau))
    return float(mesh.volumes @ P.sum(axis=1) / mesh.volumes.sum())


def _reference_pc(coord: CoordSystem, tau: float) -> float:
    P = fsg.build_cp_matrix(_NoClosure(kind=_KIND[coord]), np.array([[tau]]), np.array([1.0]),
                            np.array([_MEASURE_CONSTANT[coord]]), 1.0)
    return float(P[0, 0, 0])


#: Relative bands, 10 x `[M]` 2026-10-10 rounded up (``ta/d1/probe.py``): the sphere's production 2.1e-11 (tau = 5,
#: composite Gauss-Legendre in y on the cell edges), its reference 6.2e-13 (tau = 0.1); the slab is exact to rounding
#: (4e-16). The cylinder reads the Ki_3 Chebyshev interpolant, documented at 1e-6 absolute, which the probability
#: carries as an absolute 1.18e-6 / tau, 1.28e-6 / tau, 1.28e-6 / tau at tau = 0.1, 1, 5 (both routes): its band is
#: absolute, 4e-6 / tau, derived from the interpolant and about 3x the measurement.
_RTOL = {("production", CoordSystem.SPHERICAL): 2e-10, ("reference", CoordSystem.SPHERICAL): 1e-11,
         ("production", CoordSystem.CARTESIAN): 1e-13, ("reference", CoordSystem.CARTESIAN): 1e-13}


@pytest.mark.l1
@pytest.mark.rests_on("tests/gates/derivations/test_chord_half_lengths.py::test_each_shell_of_the_half_chord_is_exact_to_a_few_ulp",
                      "tests/gates/geometry/test_chart.py::test_the_beam_density_of_each_chart")
@pytest.mark.parametrize("route", ["production", "reference"])
@pytest.mark.parametrize("coord", [pytest.param(CoordSystem.SPHERICAL, id="sphere", marks=pytest.mark.catches("ERR-107")),
                                   pytest.param(CoordSystem.CYLINDRICAL, id="cylinder"),
                                   pytest.param(CoordSystem.CARTESIAN, id="slab")])
@pytest.mark.parametrize("tau", [0.1, 1.0, 5.0])
def test_the_first_collision_probability_is_the_closed_form(route: str, coord: CoordSystem, tau: float) -> None:
    """D1: the volume-averaged row sum of ``P_cell`` of a homogeneous unit body (production: two regions (0.5, 1.0),
    3 + 4 cells, vacuum; reference: one region) equals the closed-form ``P_c``.

    First reds (``ta/battery``): the sphere's density spelled ``y`` again (production or reference: ratio pi, every
    sphere row of that route); the slab's one-point weight back at 0.25 with the 2 gone (every reference slab row,
    ratio 2); the cylinder's density scaled by 1 + 1e-3.
    """
    got = (_production_pc if route == "production" else _reference_pc)(coord, tau)
    want = _truth(coord, tau)
    if coord is CoordSystem.CYLINDRICAL:
        band = 4e-6 / tau
        assert abs(got - want) <= band, f"{route} cylinder P_c({tau}) = {got!r}, closed form {want!r}, band {band:.1e}"
    else:
        rtol = _RTOL[(route, coord)]
        assert abs(got - want) <= rtol * want, f"{route} {_CHART[coord]} P_c({tau}) = {got!r}, closed form {want!r}"


# ── the vacuum eigenvalue (a defect left in place by ruling) ────────────────

_VACUUM_K_REASON = (
    "#546: CPSolver.compute_keff divides production by net removal with no leakage term, "
    "so solve_cp's k under a vacuum outer law is not the eigenvalue of the CP pencil (the CP twin of ERR-064); "
    "not fixed by the user's ruling of 2026-10-10 (CP's pencil is a later campaign)"
)
_MATERIALS = {2: get_mixture("A", "2g"), 0: get_mixture("B", "2g")}
#: Tight tolerances, so a k comparison reads 10 x keff_tol = 1e-10 (the white control's band).
_PARAMS = CPParams(keff_tol=1e-11, flux_tol=1e-11, max_outer=5000)


def _vacuum_mesh(coord: CoordSystem):
    """The 2G A|B body with its outer face the tag 0 (the wall it is bound to is the caller's: ``_WALLS``)."""
    geom = cp_body(coord, (0.0, 0.5, 1.0), (2, 0), 0)
    return Mesher(geom).partition((CellsByCount.uniform_volume(4), CellsByCount.uniform_volume(4))).mesh


def _pencil_k(mesh, wall) -> float:
    """The dominant eigenvalue of ``(I - T S)^{-1} T F`` from CP's own ``P_inf``: ``T`` maps a per-volume source to
    the collision-density flux, ``phi_ig = sum_j P_ji V_j q_jg / (Sigma_t,ig V_i)``; ``S`` is in-cell scattering
    (``SigS[0]`` transposed, to <- from); ``F`` is ``chi nu Sigma_f``."""
    cpm = CPMesh(mesh, boundary_conditions={0: wall})
    xs = assemble_cell_xs(_MATERIALS, mesh.mat_ids)
    n, ng = mesh.N, xs.sig_t.shape[1]
    V = mesh.volumes
    T = np.zeros((n * ng, n * ng))
    S = np.zeros_like(T)
    F = np.zeros_like(T)
    for g in range(ng):
        P = cpm.compute_pinf_group(xs.sig_t[:, g])
        T[g::ng, g::ng] = (P.T * V[None, :]) / (xs.sig_t[:, g] * V)[:, None]
    for i in range(n):
        mix = _MATERIALS[int(mesh.mat_ids[i])]
        S[i * ng:(i + 1) * ng, i * ng:(i + 1) * ng] = mix.SigS[0].toarray().T
        F[i * ng:(i + 1) * ng, i * ng:(i + 1) * ng] = np.outer(xs.chi[i], xs.sig_p[i])
    return float(np.max(np.abs(np.linalg.eigvals(np.linalg.solve(np.eye(n * ng) - T @ S, T @ F)))))


_VACUUM_PARAMS = [pytest.param(c, id=_CHART[c], marks=pytest.mark.xfail(strict=True, raises=AssertionError, reason=_VACUUM_K_REASON))
                  for c in (CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL)]


@pytest.mark.l1
@pytest.mark.rests_on(_HERE + "test_the_first_collision_probability_is_the_closed_form")
@pytest.mark.parametrize("coord", _VACUUM_PARAMS)
def test_the_vacuum_k_is_the_eigenvalue_of_the_cp_pencil(coord: CoordSystem) -> None:
    """2G A|B body (0.5, 1.0), 4 + 4 cells, vacuum outer law: ``solve_cp``'s k equals the dense pencil eigenvalue to
    10 x keff_tol relative. Strict xfail (the reason names the defect). `[M]` memo 2026-10-10: slab 1.230176 against 0.501609;
    cylinder 0.938053 against 0.041542; sphere (density fixed) about 0.59 against 0.0184. When ``compute_keff``
    carries the leakage the rows XPASS, strictness turns them red, and the marks are deleted. The marks carry
    ``raises=AssertionError`` (wave 2, S4): only the comparison failing is the expected failure, never an exception of the
    build (the sibling file ``test_cp_vacuum_not_above_white.py`` states the same defect as the physics inequality).
    """
    mesh = _vacuum_mesh(coord)
    with contextlib.redirect_stdout(io.StringIO()):
        k = solve_cp(_MATERIALS, mesh, _PARAMS, boundary_conditions={0: WallCondition.vacuum()}).keff
    reference = _pencil_k(mesh, WallCondition.vacuum())
    assert abs(k - reference) <= 10 * _PARAMS.keff_tol * reference, f"solve_cp k {k!r}, pencil eigenvalue {reference!r}"


@pytest.mark.l1
@pytest.mark.parametrize("coord", [CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL],
                         ids=["slab", "cylinder", "sphere"])
def test_under_white_the_pencil_eigenvalue_is_solve_cps_k(coord: CoordSystem) -> None:
    """The reference's own control: under a white outer law, where ``compute_keff``'s no-leakage premise holds (the
    closure returns every escaping neutron), the dense pencil eigenvalue equals ``solve_cp``'s k to 10 x keff_tol
    relative (``keff_tol = 1e-11``). This is what licenses reading the vacuum rows' disagreement as
    ``compute_keff``'s and not the pencil's.
    """
    mesh = _vacuum_mesh(coord)
    with contextlib.redirect_stdout(io.StringIO()):
        k = solve_cp(_MATERIALS, mesh, _PARAMS, boundary_conditions={0: WallCondition.white()}).keff
    reference = _pencil_k(mesh, WallCondition.white())
    assert abs(k - reference) <= 10 * _PARAMS.keff_tol * reference, f"white {_CHART[coord]}: solve_cp k {k!r}, pencil {reference!r}"


@pytest.mark.foundation
def test_the_vacuum_k_rows_are_strict_xfails() -> None:
    """The vacuum-k rows carry one ``xfail(strict=True)`` each, its reason naming the defect: a mark that lost
    ``strict`` inside ``pytest.param(marks=...)`` would turn the repair into a silent XPASS."""
    assert len(_VACUUM_PARAMS) == 3
    for param in _VACUUM_PARAMS:
        xfails = [m for m in param.marks if m.name == "xfail"]
        assert len(xfails) == 1 and xfails[0].kwargs.get("strict") is True, param.id
        assert xfails[0].kwargs.get("raises") is AssertionError, f"{param.id}: without raises=AssertionError an ImportError or a TypeError reads as the expected failure"

        assert "compute_keff" in xfails[0].kwargs.get("reason", ""), param.id
