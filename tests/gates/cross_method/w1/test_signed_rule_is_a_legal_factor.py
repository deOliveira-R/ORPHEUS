r"""G04w: a signed quadrature rule is a legal angular factor, and the vacuum-faced solve on it is unchanged (N04, spec 7.5, decision ``[ON04-0]``).

(a) The 16 tests of ``tests/gates/numerics/test_frame.py`` (11) and ``test_spherical_harmonic_space.py`` (5) the spec names exist (``ast``: each name is defined in its file); a gate cannot assert
that other tests pass, so the row is the CI step running those two files, listed with this file. (b) The frame's ``measure_space`` and ``basis_space`` on ``lebedev(13)``, ``(25)``, ``(27)`` at L = 0 and 2
build, and ``lebedev(13)`` at L = 7 is refused by the frame's door only, naming "Krein" and "#706". (c) EXACT leg: ``solve_sn`` on the vacuum-faced 2 by 2 box, one group, scattering order 0 and 1, on
``lebedev(13)`` does not raise and the solution reads converged. (d) the LOCAL value leg lives in ``test_signed_rule_keff_values.py`` (host-keyed, outside the CI step): keff on ``lebedev(13)`` equals 0.29309086 (order 0) and 0.27773609 (order 1) within the solve's own convergence bound, read from
its exit report (``Asserted(bound=1e-06, ...)``), not a chosen tolerance, and ``lebedev(11)``'s 0.29347604 and 0.27820053 for the neighbour; the values are the opening's, recomputed on the base tip by this gate's own
probe (``[M]`` 0.29347604426414825, 0.2782005320318106, 0.29309086444039906, 0.27773609108832553).

GREEN at the base tip (``[M]``): the spec says so, "on the tip all of it passes, so (a) to (c) are red only on an N04 tree that refuses a signed rule"; recorded as none at the base, with the mutation:
reinstalling the negative-weight refusal in ``Factor._require_constructible`` reddens (b) and (c) and the 16 lines of (a). Legs (a) to (c) are exact on any host (a parse, builds, a refusal and a convergence flag);
(d) compares a computation with a recorded value to a solver tolerance and is host-keyed by the spec's ruling, so it is LOCAL.
"""

from __future__ import annotations

import ast
import warnings
from pathlib import Path

import numpy as np
import pytest

from tests.gates.cross_method.w1 import _fixtures as fx

pytestmark = pytest.mark.foundation

ROOT = Path(__file__).resolve().parents[4]
FRAME_TESTS = (
    "test_analysis_hilbert_adjoint_falls_out_of_the_frame_spaces",
    "test_faces_compose_through_operator_product_with_real_spaces",
    "test_frame_faces_two_axis_contract",
    "test_measure_space_carries_the_quadrature_weights_as_its_metric",
    "test_parseval_analysis_is_an_isometry_onto_its_image",
    "test_parseval_frame_square_closes",
    "test_petrov_galerkin_degenerate_equals_galerkin_bit_identically",
    "test_petrov_galerkin_degenerate_project_equals_galerkin_project",
    "test_pi_R_is_4pi_identity_through_the_frame",
    "test_reconstruction_hilbert_adjoint_falls_out_of_the_frame_spaces",
    "test_the_frame_can_now_spell_its_own_projector",
)
SH_TESTS = (
    "test_H_equals_parseval_metric_times_S0",
    "test_R_transpose_carries_d_ell_and_RH_carries_d_ell_squared",
    "test_R_transpose_equals_2l_plus_1_times_S0_transpose",
    "test_T_carries_w_n_and_H_carries_the_parseval_metric",
    "test_pi_R_is_4pi_identity_on_band_limited",
)


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_the_sixteen_named_frame_lines_exist() -> None:
    for rel, names in (("tests/gates/numerics/test_frame.py", FRAME_TESTS), ("tests/gates/numerics/test_spherical_harmonic_space.py", SH_TESTS)):
        defined = {n.name for n in ast.walk(ast.parse((ROOT / rel).read_text())) if isinstance(n, ast.FunctionDef)}
        _require(len(names) == (11 if "test_frame" in rel else 5), f"{rel}: the named list is not the spec's")
        _require(set(names) <= defined, f"{rel}: named tests that are not defined: {sorted(set(names) - defined)}")


def test_the_frame_spaces_build_on_a_signed_rule_and_the_door_refuses_the_under_resolved_basis() -> None:
    from orpheus.numerics.frame_faces import angular_frame, basis_space, measure_space
    from orpheus.numerics.quadrature.directional import Quadrature

    for order in (13, 25, 27):
        for L in (0, 2):
            frame = angular_frame(Quadrature.lebedev(order), L)
            _require(measure_space(frame) is not None and basis_space(frame) is not None, f"lebedev({order}) L={L}: a space did not build")
    try:
        basis_space(angular_frame(Quadrature.lebedev(13), 7))
    except ValueError as e:
        _require("Krein" in str(e) and "#706" in str(e), f"lebedev(13) L=7 is refused, but not by the frame's door naming Krein and #706: {e}")
    else:
        pytest.fail("lebedev(13) at L = 7 was accepted")


def _vacuum_box_solve(order: int, scattering_order: int):
    from orpheus.mesh import Mesh2D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn
    from orpheus.transport.boundary_conditions import WallCondition

    m0, m1, _ = fx._sn_mixtures(1)
    mat_map = np.zeros((2, 2), dtype=int)
    mat_map[:1, :] = 2
    mesh = Mesh2D(edges_x=np.linspace(0, 1.3, 3), edges_y=np.linspace(0, 0.7, 3), mat_map=mat_map, face_roles={"xmin": 0, "xmax": 0, "ymin": 0, "ymax": 0})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return solve_sn(
            {2: m0, 0: m1},
            mesh,
            Quadrature.lebedev(order),
            scattering_order=scattering_order,
            boundary_conditions=fx.conditions_for(mesh, {0: WallCondition.vacuum()}),
            max_outer=300,
        )


@pytest.mark.parametrize("scattering_order", [0, 1])
def test_the_vacuum_faced_solve_on_lebedev_13_converges(scattering_order: int) -> None:
    solution = _vacuum_box_solve(13, scattering_order)
    _require(bool(solution.converged), f"the solve on lebedev(13) at scattering order {scattering_order} does not read converged")
    _require(np.isfinite(float(solution.outcome.keff)) and float(solution.outcome.keff) > 0.0, "keff is not a positive finite number")
