r"""G03u: a reflective solve on a signed rule is refused upstream, so no raise is needed at the kernel projection gauge (N03b, spec 7.4b).

A 2 by 2 box, every face specular with albedo 1, one group, two materials: ``solve_sn`` on ``Quadrature.lebedev(13)`` (a rule with 8 negative weights)
raises ``ValueError`` whose message carries the fragment "carry a vanishing metric", and ``solve_sn`` on ``lebedev(11)`` solves. The refusal comes
from ``OrthonormalSpanBasis.of_generators`` (the gate reads the traceback's frames), upstream of ``KernelProjectionGauge.displacement``, which is
why N03b adds no raise there: a guard nothing reaches is not written. The message names a vanishing metric for what is a negative weight; the
wording and the signed orthonormalisation are open under #706.

Green at the base of N03b (the row pins a refusal that exists). Mutation: ``OrthonormalSpanBasis.of_generators`` with its weight check removed (a
copy of the function that skips the ``metric > 0`` refusal, bound in every module that holds the class) must redden the gate: the solve on
``lebedev(13)`` then raises something else or solves. Exact on any host (a refusal, not a value).
"""

from __future__ import annotations

import traceback
import warnings
from typing import Any

import numpy as np
import pytest

from tests.gates.cross_method.w1 import _fixtures as fx

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _reflective_solve(order: int) -> Any:
    from orpheus.mesh import Mesh2D
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn
    from orpheus.transport.boundary_conditions import WallCondition

    m0, m1, _ = fx._sn_mixtures(1)
    mat_map = np.zeros((2, 2), dtype=int)
    mat_map[:1, :] = 2
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1.3, 3),
        edges_y=np.linspace(0, 0.7, 3),
        mat_map=mat_map,
        face_roles={"xmin": 0, "xmax": 0, "ymin": 0, "ymax": 0},
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return solve_sn(
            {2: m0, 0: m1},
            mesh,
            Quadrature.lebedev(order),
            scattering_order=0,
            boundary_conditions=fx.conditions_for(mesh, {0: WallCondition.specular(1.0)}),
            max_outer=300,
        )


def test_a_reflective_solve_on_a_signed_rule_is_refused_upstream() -> None:
    solution = _reflective_solve(11)
    keff = float(solution.outcome.keff)
    _require(np.isfinite(keff) and keff > 0.0, f"the reflective solve on lebedev(11) reads keff {keff}")
    try:
        _reflective_solve(13)
    except ValueError as e:
        _require("carry a vanishing metric" in str(e), f"the refusal does not carry the fragment 'carry a vanishing metric': {e}")
        frames = traceback.extract_tb(e.__traceback__)
        origin = [f for f in frames if f.name == "of_generators" and f.filename.endswith("orthonormal_span_basis.py")]
        _require(
            len(origin) == 1,
            f"the refusal does not originate in OrthonormalSpanBasis.of_generators: frames {[(f.name, f.filename.rsplit('/', 1)[-1]) for f in frames[-4:]]}",
        )
        _require(
            not any(f.name == "displacement" for f in frames),
            "the refusal passed through KernelProjectionGauge.displacement, which this unit leaves unguarded",
        )
    else:
        pytest.fail("a reflective solve on lebedev(13) was not refused")
