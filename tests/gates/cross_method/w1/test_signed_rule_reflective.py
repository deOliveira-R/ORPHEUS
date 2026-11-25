r"""G03u: a reflective solve on a signed rule is refused upstream, so no raise is needed at the kernel projection gauge (N03b, spec 7.4b).

A 2 by 2 box, every face specular with albedo 1, one group, two materials: ``solve_sn`` on ``Quadrature.lebedev(13)`` (a rule with 8 negative weights)
raises ``ValueError`` whose message carries the fragments "Krein" and "#706" (ruling Q1 of the orchestrator, 2026-11-25), and ``solve_sn`` on ``lebedev(11)`` solves. The refusal comes
from ``OrthonormalSpanBasis.of_generators`` (the gate reads the traceback's frames), upstream of ``KernelProjectionGauge.displacement``, which is
why N03b adds no raise there: a guard nothing reaches is not written. The refusal is the one a verb that needs a positive form
makes (``_require_positive``, the helper of ``FunctionSpace.norm``), through the one predicate ``HilbertMetric.positive``; a ZERO weight keeps its
own refusal, "carry a vanishing metric", which is true of a zero (no minimum-norm member, a division by zero), and a second row pins it.

Green at the base of N03b (``c18f156f``): the row pins a refusal that already exists. Mutations: the negative case routed back to the old message reddens the
fragments; ``HilbertMetric.positive`` forced true removes the refusal and reddens the gate; the zero case dropped reddens the zero-weight row.
The original mutation, ``OrthonormalSpanBasis.of_generators`` with its refusal removed, reddens the gate: the solve on ``lebedev(13)`` then raises something else (the frame's own Gram door, behind it) or
solves. Exact on any host (a refusal, not a value).
"""

from __future__ import annotations

import re
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
        _require("Krein" in str(e) and "#706" in str(e), f"the refusal does not carry the fragments 'Krein' and '#706': {e}")
        # the count is of the SUPPORT ROWS of the refused summand that carry a negative metric (ordinates times face cells), so it
        # depends on the mesh; the gate asks only that the message states a positive count
        counted = re.search(r"has (\d+) negative direction", str(e))
        _require(counted is not None and int(counted.group(1)) > 0, f"the refusal does not state a positive count of negative support rows: {e}")
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


def test_a_zero_weight_keeps_its_own_refusal_and_a_negative_one_is_refused_as_a_verb_needing_a_positive_form() -> None:
    """The two cases the refusal conflated, called directly: a zero weight is a vanishing metric (no minimum-norm member), a negative weight needs a Krein arm."""
    from orpheus.numerics.basis.orthonormal_span_basis import OrthonormalSpanBasis

    columns = np.ones((3, 1))
    try:
        OrthonormalSpanBasis.of_generators(columns, np.array([1.0, 0.0, 2.0]), label="zero")
    except ValueError as e:
        _require("carry a vanishing metric" in str(e), f"a zero weight does not carry the fragment 'carry a vanishing metric': {e}")
        _require("Krein" not in str(e), f"a zero weight is refused as a Krein case: {e}")
    else:
        pytest.fail("a zero weight was not refused")
    try:
        OrthonormalSpanBasis.of_generators(columns, np.array([1.0, -1.0, 2.0]), label="negative")
    except ValueError as e:
        _require("Krein" in str(e) and "#706" in str(e) and "1 negative direction" in str(e), f"a negative weight does not name the Krein arm, #706 and the count: {e}")
        _require("vanishing" not in str(e), f"a negative weight is refused with the false reason of a vanishing metric: {e}")
    else:
        pytest.fail("a negative weight was not refused")
