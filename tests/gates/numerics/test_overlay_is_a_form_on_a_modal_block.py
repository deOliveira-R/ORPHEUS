r"""A dense form lives only on a MODAL axis in production, and no form can sit beside the counting energy axis.

The space-membership unit U2 (plan ``.claude/plans/space_membership_covector_gauge.md``) refused a positioned dense form on a nodal axis, which
carries the counting measure, and on the energy axis, where it had read 83 in a space's pairing against the reaction rate's 46. Since unit N04
the form is the factor's own and no overlay exists to refuse: an energy axis has the form ``NoForm``, so a factor on it takes a
``FormlessVariance`` and a formed variance is refused by type (the leg below), and a coefficient (Modal) axis holding a dense form is a factor on that form,
which the space pairs by (the admitted leg and the production dressing).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis, Modal
from orpheus.numerics.declaration import Coord
from orpheus.numerics.factor import Factor, FormlessVariance, Kind, Subquotient, Variance, VarianceChart
from orpheus.numerics.metric import DenseMetric, Face
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import TensorSpace

pytestmark = pytest.mark.l0

_GRAM = DenseMetric(np.array([[2.0, 0.5], [0.5, 1.0]]), kernel=Coord(2, ()), face=Face.METRIC)


@pytest.mark.catches("ERR-121")
def test_a_form_on_the_counting_energy_axis_is_REFUSED():
    """The energy-block leak, closed by type: the energy axis's form is ``NoForm``, so a factor on it with a formed variance is a ``TypeError``,
    and with its ``FormlessVariance`` it holds no form for the walk to apply."""
    energy = EnergyAxis.synthetic(2)
    with pytest.raises(TypeError, match="formless"):
        Factor(energy, Variance(Kind.MEASURE, VarianceChart.DENSITY), Subquotient.whole(2))
    space = TensorSpace.of_axes(Factor(energy, FormlessVariance.MEASURE, Subquotient.whole(2)))
    x = np.array([3.0, -1.0])
    assert space.apply_metric(x).tobytes() == x.tobytes()


def test_a_dense_form_on_a_declared_nodal_axis_is_that_axis_form():
    """A declared nodal axis holding a dense form is a factor on that form (the first red of the retired refusal is now the unit's own admitted case)."""
    nodal = Axis("a", (2,), origin=Declared(BasisKind.NODAL, _GRAM))
    space = TensorSpace.of_axes(Factor(nodal, Variance(Kind.MEASURE, VarianceChart.DENSITY), Subquotient.whole(2)))
    x, y = np.array([1.0, 2.0]), np.array([3.0, -1.0])
    by_hand = float(x @ np.array([[2.0, 0.5], [0.5, 1.0]]) @ y)
    assert space.inner_product(x, y) == pytest.approx(by_hand, rel=4 * np.finfo(float).eps)


def test_a_form_on_a_modal_axis_is_the_metric():
    """The coefficient axis: the form IS the block's metric, read by the pairing."""
    coefficients = Axis("a", (2,), origin=Declared(Modal(np.array([1.0, 0.0])), _GRAM))
    space = TensorSpace.of_axes(Factor(coefficients, Variance(Kind.MEASURE, VarianceChart.DENSITY), Subquotient.whole(2)))
    x, y = np.array([1.0, 2.0]), np.array([3.0, -1.0])
    by_hand = float(x @ np.array([[2.0, 0.5], [0.5, 1.0]]) @ y)
    assert space.inner_product(x, y) == pytest.approx(by_hand, rel=4 * np.finfo(float).eps)


def test_the_production_dense_dressing_is_admitted():
    """Positive control: the production mint, ``basis_space``, makes the head ``FromFrame``, whose kind is Modal, and the dense Parseval form of a
    DENSE-verdict frame is the form of its factor."""
    head = frame_faces.basis_space(frame_faces.angular_frame(Quadrature.folded_product(4, 4), 2))
    assert isinstance(head.axes[0].kind, Modal)
    assert isinstance(head.factors[0].form, DenseMetric), "premise: this frame's dressing is the dense arm"
