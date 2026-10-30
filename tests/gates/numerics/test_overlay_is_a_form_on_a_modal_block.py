r"""A positioned form lives only on a MODAL axis: a nodal axis carries a measure, the counting measure included.

The space-membership unit U2 (plan ``.claude/plans/space_membership_covector_gauge.md``, spec
``scratch/membership_covector/spec/spec.md`` §U2.1). A dense form is a Gram on expansion COEFFICIENTS, so it is
positioned only on an axis whose kind is :class:`~orpheus.numerics.axis.Modal`. A NODAL axis with ``weights=None`` is
not measure-free: it carries the counting measure, and a form beside it is a second metric source on one block.

Before U2 the door (``_require_positioned_over``) refused a form only on an axis that carried explicit weights, so a
form on a counting axis was admitted. `[M]` 2026-10-30, the spec's overlay census (31 test files that touch dense
metrics, 1264 passed, every admitted positioned form recorded by its axis): 137 forms admitted, 108 of them on the
production frame-dressed heads (``FromFrame``, Modal) and 0 of the 22 the predicate reclassifies in production. The
one latent defect it closes is the energy-block overlay: a dense form on the counting energy axis of a cross-section
space read 83 in the space's pairing where the reaction rate, which counts groups, read 46.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics import frame_faces
from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis, Modal
from orpheus.numerics.metric import DenseMetric, FactoredMetric
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.l0

_E = EnergyAxis.synthetic(2)
_S = Axis("spatial", (3,), weights=np.array([1.0, 2.0, 3.0]), origin=Declared(BasisKind.NODAL))
_GRAM = DenseMetric(np.array([[2.0, 0.5], [0.5, 1.0]]))


def test_a_form_on_the_counting_energy_axis_is_REFUSED():
    """The energy-block overlay leak. First red: admitted before U2, and the pairing then read 83 against the reaction
    rate's 46 on one flux."""
    with pytest.raises(ValueError, match="counting measure"):
        FunctionSpace("xs", (2, 3), axes=(_E, _S), metric=FactoredMetric((((2,), _GRAM), ((3,), None))))


def test_a_form_on_a_counting_nodal_axis_is_REFUSED():
    """The same law on a declared nodal axis. First red: ``test_space_of_axes.py``'s B2 row admitted it."""
    counting = Axis("a", (2,), origin=Declared(BasisKind.NODAL))
    with pytest.raises(ValueError, match="counting measure"):
        FunctionSpace("bad", (2,), axes=(counting,), metric=FactoredMetric((((2,), _GRAM),)))


def test_a_form_on_a_modal_axis_is_admitted_and_is_the_metric():
    """The admitted leg: on a coefficient axis the form IS the block's metric, read by the pairing."""
    coefficients = Axis("a", (2,), origin=Declared(Modal(np.array([1.0, 0.0]))))
    space = FunctionSpace("good", (2,), axes=(coefficients,), metric=FactoredMetric((((2,), _GRAM),)))
    x, y = np.array([1.0, 2.0]), np.array([3.0, -1.0])
    by_hand = float(x @ np.array([[2.0, 0.5], [0.5, 1.0]]) @ y)
    assert space.inner_product(x, y) == pytest.approx(by_hand, rel=4 * np.finfo(float).eps)


def test_the_production_dense_dressing_is_admitted():
    """Positive control: the production mint, ``frame_faces._dressed`` through ``basis_space``, makes the head
    ``FromFrame``, whose kind is Modal, so the dense Parseval form of a DENSE-verdict frame is admitted. A predicate
    that refused every form would redden this row."""
    head = frame_faces.basis_space(frame_faces.angular_frame(Quadrature.folded_product(4, 4), 2))
    assert head.axes is not None and isinstance(head.axes[0].kind, Modal)
    assert isinstance(head.metric, FactoredMetric), "premise: this frame's dressing is the dense arm"
