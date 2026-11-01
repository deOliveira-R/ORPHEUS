"""MOVED from tests/gates/numerics/test_angular_symmetry.py (layer DAG Phase 2 refresh, specs_2r §5). Bodies unchanged."""

from __future__ import annotations
import pathlib
import numpy as np
import pytest
from scipy.sparse import csr_matrix
from orpheus.data.macro_xs.mixture import Mixture
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.quadrature.registry import GEOMETRY_ANGULAR_SYMMETRY
from orpheus.geometry.symmetry import SubgroupOfO3
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for


_CONDITIONS = {0: WallCondition.vacuum()}


ORPHEUS_ROOT = pathlib.Path(__file__).resolve().parents[3]


S = SubgroupOfO3
































def _one_group_mixture(sigma_t: float = 1.0, sigma_s: float = 0.5) -> Mixture:
    zero = np.zeros(1)
    return Mixture(
        SigC=np.array([sigma_t - sigma_s]), SigL=zero.copy(), SigF=zero.copy(),
        SigP=zero.copy(), SigT=np.array([sigma_t]),
        SigS=[csr_matrix(np.array([[sigma_s]]))],
        Sig2=[csr_matrix(np.zeros((1, 1)))], chi=zero.copy(),
    )


def _ordinate_permutation(quad: Quadrature, matrix: np.ndarray, atol: float = 1e-10):
    """``n -> n'`` with ``Omega_{n'} = M Omega_n``, or ``None`` if the rule is not
    closed under ``M``."""
    omega = np.column_stack([quad.mu_x, quad.mu_y, quad.mu_z])
    mapped = omega @ np.asarray(matrix).T
    out = []
    for row in mapped:
        distance = np.abs(omega - row).max(axis=1)
        j = int(np.argmin(distance))
        if distance[j] > atol:
            return None
        out.append(j)
    perm = np.asarray(out)
    return perm if len(set(out)) == len(out) else None


def _interior(solution) -> np.ndarray:
    field = solution.angular_flux
    return np.asarray(field.interior.values if hasattr(field, "interior") else field.values)


def _evenness(psi: np.ndarray, perm: np.ndarray) -> float:
    axis = 0 if psi.shape[0] == len(perm) else 1
    return float(np.abs(np.take(psi, perm, axis=axis) - psi).max()
                 / max(np.abs(psi).max(), 1e-300))


@pytest.fixture(scope="module")
def _cylinder_solution():
    quad = Quadrature.folded_product(4, 8)
    mesh = Mesher(StructuredGeometry.cylinder((0.0, 2.0), (0,), outer=0),
                  ).partition(CellsByCount.uniform_width(8)).mesh
    source = np.zeros((quad.n_ordinates, 1, 8))
    source[:, :, :4] = 1.0                     # r-ASYMMETRIC: no flat-flux nulling
    solution = solve_sn_fixed_source({0: _one_group_mixture()}, mesh, quad, source,
                                     inner_tol=1e-13, max_inner=20_000, boundary_conditions=conditions_for(mesh, _CONDITIONS))
    return quad, _interior(solution)


@pytest.fixture(scope="module")
def _plane_solution():
    quad = Quadrature.product(4, 8)
    mesh = Mesh2D(edges_x=np.linspace(0.0, 2.0, 6), edges_y=np.linspace(0.0, 1.0, 4),
                  mat_map=np.zeros((5, 3), dtype=int), coord=CoordSystem.CARTESIAN,
                  face_roles={
                      "xmin": 0,
                      "xmax": 0,
                      "ymin": 0,
                      "ymax": 0,
                  })
    source = np.zeros((quad.n_ordinates, 1, 5, 3))
    source[:, :, :2, :1] = 1.0                 # x- AND y-asymmetric corner source
    solution = solve_sn_fixed_source({0: _one_group_mixture()}, mesh, quad, source,
                                     inner_tol=1e-13, max_inner=20_000, boundary_conditions=conditions_for(mesh, _CONDITIONS))
    return quad, _interior(solution)


@pytest.fixture(scope="module")
def _slab_solution():
    quad = Quadrature.gauss_legendre(8)
    mesh = Mesher(StructuredGeometry.slab((0.0, 2.0), (0,), left=0,
                                          right=0),
                  ).partition(CellsByCount.uniform_width(8)).mesh
    source = np.zeros((quad.n_ordinates, 1, 8))
    source[:, :, :4] = 1.0
    solution = solve_sn_fixed_source({0: _one_group_mixture()}, mesh, quad, source,
                                     inner_tol=1e-13, max_inner=20_000, boundary_conditions=conditions_for(mesh, _CONDITIONS))
    return quad, _interior(solution)


_MIRRORS = {"sigma_x": np.diag([-1.0, 1.0, 1.0]),
            "sigma_y": np.diag([1.0, -1.0, 1.0]),
            "sigma_z": np.diag([1.0, 1.0, -1.0])}


class TestR3UnspentIsMeasuredOnTheSolver:
    r"""**C1-C3** — ``unspent`` is a claim about the SOLUTION, and this is the most a
    gate can say about it.

    ⛔ **What a gate CANNOT do.**  ``unspent`` is a statement about the CONTINUUM
    transport equation of a geometry — "psi is even under this finite group acting on
    directions alone".  No test verifies that; it is a derivation (the archivist's
    page), and a wrong derivation would make the geometry's whole reduction wrong,
    not just its registry row.  A test cannot even enumerate the claim's domain: it
    ranges over every source, every material, every boundary law.

    ⭐ **What a gate CAN do, and it is the load-bearing half.**  Every downstream fold
    rests on the DISCRETE solver's answer, and there the claim is exactly checkable:
    solve a deliberately asymmetric fixed-source problem and compare ``psi`` at
    ordinate ``n`` with ``psi`` at the ordinate the group element maps it to.  A
    positive leg per element of ``unspent``, a NEGATIVE leg per element outside it —
    without the second the row is `vv-principles` #19-blind, because a solver whose
    psi is even under EVERYTHING would pass the positive leg too.

    `[M]` 2026-09-03, all three fixtures in **1.0 s** total:

    ==============  ================  =====================  ==================
    geometry        rule              ``sigma_z``            ``sigma_y``
    ==============  ================  =====================  ==================
    cartesian2d     product(4,8)      **0.0** (EVEN)         6.043e-01 (NOT)
    cylinder        folded_product    **0.0** (EVEN)         no permutation
    slab            gauss_legendre    0.0 but VACUOUS        0.0 but VACUOUS
    ==============  ================  =====================  ==================

    and ``sigma_x`` moves every one of them: 8.175e-01 / 7.224e-01 / 6.493e-01.

    ⚠ **The trap the vacuity guard caught.**  On the SLAB the ``sigma_y`` and
    ``sigma_z`` permutations are the IDENTITY (a 1-D polar rule carries
    ``mu_y = mu_z = 0`` for every ordinate — ERR-080's territory), so their
    ``0.0`` readings are `vv-principles` Mode-8 tautologies, not evenness.  Reported
    without the guard they would have contradicted ``unspent = Trivial``.  Every
    positive leg here therefore asserts the permutation MOVES ordinates.
    """

    @pytest.mark.foundation
    def test_c1_cartesian2d_is_even_under_sigma_z_and_under_NOTHING_else(
        self, _plane_solution
    ) -> None:
        quad, psi = _plane_solution
        plane = GEOMETRY_ANGULAR_SYMMETRY["cartesian2d"]
        assert plane.unspent == S.Mirror("z")

        perm = _ordinate_permutation(quad, _MIRRORS["sigma_z"])
        assert perm is not None and int(np.count_nonzero(perm != np.arange(len(perm)))) == 32
        assert _evenness(psi, perm) == 0.0, "z-uniformity is EXACT in the 2-D operator"

        for name in ("sigma_x", "sigma_y"):
            other = _ordinate_permutation(quad, _MIRRORS[name])
            assert other is not None
            assert int(np.count_nonzero(other != np.arange(len(other)))) == 24
            assert _evenness(psi, other) > 0.5, (
                f"{name} must MOVE the solution — it is outside cartesian2d's unspent "
                f"group, and a fold by it empties two sweep quadrants (D1)"
            )

    @pytest.mark.foundation
    def test_c2_the_cylinder_is_even_under_sigma_z_and_sigma_x_moves_it(
        self, _cylinder_solution
    ) -> None:
        r"""⛔ **The ``sigma_y`` half of ``D_1h`` is UN-WITNESSABLE on the solver
        path, by construction.**  `[M]` 2026-09-03 a cylindrical ``SNProblem`` admits
        only CARRYING quadratures (``assert_carrying_quadrature``): **15 of 15**
        ``folded_product`` rules pass, **0 of 20** ``product`` / ``lebedev`` /
        ``level_symmetric`` rules do.  Every admissible rule IS the ``sigma_y``
        quotient, so no cylindrical solve ever stores both signs of ``mu_y`` and no
        test can ask whether psi is even under the mirror the rule already spent.

        ⟹ ``sigma_y in unspent`` is, for the cylinder, a STRUCTURAL commitment the
        tree already makes one layer down — the shipped solver refuses every
        unfolded rule — and the assertion below is the strongest available form of
        it.  Read with `vv-principles` Mode 12: the functional (a folded solve)
        annihilates the error class, so no fixture can be found; saying so is the
        honest deliverable, not a green row."""
        quad, psi = _cylinder_solution
        cylinder = GEOMETRY_ANGULAR_SYMMETRY["cylinder"]
        assert cylinder.unspent == S.Dnh(1)
        assert cylinder.unspent.contains(S.Mirror("y"))
        assert cylinder.unspent.contains(S.Mirror("z"))
        assert not cylinder.unspent.contains(S.Mirror("x"))

        # the sigma_y half: no permutation exists — the rule already spent it
        assert _ordinate_permutation(quad, _MIRRORS["sigma_y"]) is None
        assert quad.measure.quotient_group == S.Mirror("y")

        # the sigma_z half: measurable, and EXACT
        perm = _ordinate_permutation(quad, _MIRRORS["sigma_z"])
        assert perm is not None and int(np.count_nonzero(perm != np.arange(len(perm)))) == 16
        assert _evenness(psi, perm) == 0.0

        # the negative leg: sigma_x is outside D_1h and moves the solution 72 %
        other = _ordinate_permutation(quad, _MIRRORS["sigma_x"])
        assert other is not None
        assert _evenness(psi, other) > 0.5

    @pytest.mark.foundation
    def test_c3_the_slab_keeps_NOTHING_unspent(self, _slab_solution) -> None:
        r"""``unspent = Trivial`` says a slab solution is even under no non-trivial
        mirror.  ``sigma_x`` (the OWED closure) moves it by 65 %, and the other two
        mirrors act as the IDENTITY on a polar rule — so their zero readings carry
        no information and are asserted VACUOUS rather than read as evenness."""
        quad, psi = _slab_solution
        slab = GEOMETRY_ANGULAR_SYMMETRY["slab"]
        assert slab.unspent == S.Trivial
        assert slab.owed == S.Mirror("x")

        perm = _ordinate_permutation(quad, _MIRRORS["sigma_x"])
        assert perm is not None
        assert int(np.count_nonzero(perm != np.arange(len(perm)))) == 8
        assert _evenness(psi, perm) > 0.5, (
            "a slab solution is not even in mu for a general source — which is why "
            "nothing finite is left UNSPENT and sigma_x is only OWED"
        )
        for name in ("sigma_y", "sigma_z"):
            vacuous = _ordinate_permutation(quad, _MIRRORS[name])
            assert vacuous is not None
            assert int(np.count_nonzero(vacuous != np.arange(len(vacuous)))) == 0, (
                f"{name} now MOVES ordinates on a polar rule — its 0.0 evenness "
                f"reading has stopped being a tautology and needs a real gate"
            )


