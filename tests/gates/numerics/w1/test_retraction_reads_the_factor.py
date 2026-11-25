r"""G04m: the retraction over an axis reads the factor's axis, restricted where the space was built restricted (N04, spec 7.5).

On a ``TensorSpace`` whose angular axis is ``Axis.restricted((0, 2))`` the DOMAIN of ``AxisRetractionOperator.of(s, "angular")`` has angular extent 2 (shape ``(2, 2, 6)``) and the
codomain is the marginal ``(2, 6)``.

First red at the base tip: import-only (``TensorSpace`` absent; ``FunctionSpace.of_axes`` already reads the restricted axis, ``(2, 2, 6)`` on the tip). Mutation: building the
space on the UNRESTRICTED axis in the gate reddens the domain shape (the parent-axis mutation is unspellable: ``Restriction`` holds the parent's digest and the index set, never
the parent axis). Exact on any host.
"""

from __future__ import annotations

import pytest

from tests.gates.numerics.w1 import _tensor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_retraction_reads_the_factor() -> None:
    from orpheus.numerics.axis_marginal import AxisRetractionOperator

    TS = P.tensor_space()
    angular = P.diag_axis(P.ANGULAR_W, "angular").restricted((0, 2))
    s = TS.of_axes(P.density_factor(angular), P.formless_factor(P.formless_axis(2, "energy")), P.density_factor(P.diag_axis((1.0, 2.0, 3.0, 4.0, 5.0, 6.0), "spatial")))
    op = AxisRetractionOperator.of(s, "angular")
    _require(tuple(op.domain.shape) == (2, 2, 6), f"the domain reads {tuple(op.domain.shape)}, not (2, 2, 6)")
    _require(tuple(op.codomain.shape) == (2, 6), f"the codomain reads {tuple(op.codomain.shape)}, not (2, 6)")
