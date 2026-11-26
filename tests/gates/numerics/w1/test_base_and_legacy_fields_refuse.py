r"""G04x: the abstract base does not construct, a tensor space answers no legacy metric field, a leaf states one metric source, a mixed product keeps the default carrier (N04, rulings W1, W2, W5a, W6).

(i) ``FunctionSpace`` itself is refused with a ``TypeError`` (a real ``raise``: this file runs under ``-O``) whose message says a space is a ``TensorSpace`` of factors or one of the
listed arms; each class of the allow-list ``FUNCTION_SPACE_ARMS`` still constructs, from the fixtures of ``tests/gates/test_content_identity_census.py`` (one instance per arm; the
common base ``_RadialCharacteristicSubSpace`` is reached through its two subclasses). Mutation: the refusal removed from ``FunctionSpace.__post_init__`` reddens the bare row.
(ii) on a ``TensorSpace`` the two leaf fields ``inner_product_weights`` and ``metric`` raise ``AttributeError`` naming ``factor(label).form``, and on a leaf arm they answer. Mutation: the
fields returning ``None`` on a tensor space again reddens the refusal row.
(iii) a leaf arm built with two metric sources (the array and the object) is refused, "one metric source only". Mutation: the refusal removed reddens the row.
(iv) the product of two spaces whose carriers disagree (a dual times a primal: a mixed tensor) does not raise and takes the default carrier, so it mints ``Member``; the product of two spaces
whose carriers agree keeps that carrier. Mutation: the product taking the left factor's carrier reddens the mixed row.

Exact on any host: classes, identities and a raise.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.numerics.w1 import _tensor_probe as T
from orpheus.numerics.direct_sum import BULK, DirectSumSpace, TRACE

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def arm_instances() -> dict[str, object]:
    r"""One instance per allow-listed arm, keyed by the arm's module-qualified name; the radial base is represented by its two subclasses."""
    from orpheus.numerics.axis import EnergyAxis
    from orpheus.numerics.coproduct_layout import CoproductLayout
    from orpheus.numerics.factor import Factor, Kind, VarianceChart
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.numerics.space import TensorSpace
    from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
    from orpheus.numerics.direct_sum import DirectSumSpace
    from orpheus.numerics.spaces.radial_characteristic_space import RadialCharacteristicBoundarySpace, RadialCharacteristicInteriorSpace
    from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace
    from orpheus.numerics.direct_sum import DirectSumSpace

    quad = Quadrature.gauss_legendre(4)
    layout = CoproductLayout.from_named_shapes([("xmin", (quad.n_ordinates, 1)), ("xmax", (quad.n_ordinates, 1))])
    angular = AngularTraceSpace.from_quadrature_and_layout(quad, layout)
    scalar = ScalarTraceSpace.for_faces([("xmin", ()), ("xmax", ())], 2, {"xmin": 1.0, "xmax": 1.0})
    levels, volumes = (0, 2, 5), np.array([1.5, 2.5, 3.5, 4.5])
    bulk = TensorSpace.of_axes(Factor.whole_on(EnergyAxis.synthetic(2), Kind.MEASURE, VarianceChart.DENSITY))
    instances: list[object] = [
        bulk,
        angular,
        angular.outflow_space("xmin"),
        scalar,
        RadialCharacteristicInteriorSpace.for_levels(levels, ng=2, nx=4, cell_volumes=volumes),
        RadialCharacteristicBoundarySpace.for_levels(levels, ng=2, nx=4, cell_volumes=volumes),
        DirectSumSpace.of_pairs(((BULK, bulk), (TRACE, scalar))),
        DirectSumSpace.of_pairs(tuple((f's{i}', m) for i, m in enumerate((bulk, bulk.dual())))),
        angular.dual(),
    ]
    return {f"{type(x).__module__}.{type(x).__qualname__}": x for x in instances}


def arms_missing(instances: dict[str, object], listed: object) -> list[str]:
    r"""The allow-listed classes that no instance represents (an instance of a class, or of a subclass of it, represents it)."""
    classes = [type(x) for x in instances.values()]
    missing = []
    for key in listed:  # type: ignore[attr-defined]
        represented = any(f"{b.__module__}.{b.__qualname__}" == key for c in classes for b in c.__mro__)
        if not represented:
            missing.append(key)
    return missing


def test_the_base_refuses_direct_construction() -> None:
    from orpheus.numerics.space import FUNCTION_SPACE_ARMS, FunctionSpace

    with pytest.raises(TypeError, match="not constructed directly"):
        FunctionSpace("x", (3,), inner_product_weights=np.ones(3))  # type: ignore[abstract]
    with pytest.raises(TypeError, match="TensorSpace of factors"):
        FunctionSpace("y", (2,))  # type: ignore[abstract]
    instances = arm_instances()
    _require(all(isinstance(x, FunctionSpace) for x in instances.values()), "an arm instance is not a FunctionSpace")
    _require(all(type(x) is not FunctionSpace for x in instances.values()), "an arm instance is the bare base")
    missing = arms_missing(instances, FUNCTION_SPACE_ARMS)
    _require(not missing, f"{len(missing)} of {len(FUNCTION_SPACE_ARMS)} allow-listed arms have no constructing fixture: {missing}")
    _require(arms_missing({}, ["a.B"]) == ["a.B"], "the census does not find an arm with no instance")


def test_a_tensor_space_answers_no_legacy_metric_field() -> None:
    space = T.triple_space()
    for name in ("inner_product_weights", "metric"):
        with pytest.raises(AttributeError, match=r"factor\(label\)\.form"):
            getattr(space, name)
        with pytest.raises(AttributeError):
            getattr(space.dual(), name)
    instances = arm_instances()
    angular = next(x for k, x in instances.items() if k.endswith(".AngularTraceSpace"))
    dual = next(x for k, x in instances.items() if k.endswith(".DualSpace"))
    _require(np.asarray(angular.inner_product_weights).shape == angular.shape, "a leaf arm does not answer inner_product_weights")  # type: ignore[attr-defined]
    _require(angular.metric is None and dual.metric is None and dual.inner_product_weights is None, "a leaf arm does not answer its metric fields")  # type: ignore[attr-defined]


def test_a_leaf_states_one_metric_source() -> None:
    from orpheus.numerics.metric import DiagonalMetric
    from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace

    scalar = ScalarTraceSpace.for_faces([("xmin", ()), ("xmax", ())], 2, {"xmin": 1.0, "xmax": 1.0})
    weights = np.ones(scalar.shape)
    with pytest.raises(ValueError, match="one metric source only"):
        ScalarTraceSpace(name="two", shape=scalar.shape, inner_product_weights=weights, metric=DiagonalMetric(weights), layout=scalar.layout)
    one = ScalarTraceSpace(name="one", shape=scalar.shape, inner_product_weights=weights, layout=scalar.layout)
    _require(one.metric is None, "the control with one source does not construct")


def test_a_mixed_product_takes_the_default_carrier() -> None:
    from orpheus.numerics.space import Covector, Member

    primal = T.triple_space().sub("angular")
    other = T.triple_space().sub("spatial")
    column = np.ones(primal.shape + other.shape)
    mixed = primal.dual() * other
    _require(mixed.carrier is None, "a mixed product does not take the default carrier")
    _require(type(mixed.member(column)) is Member, "a mixed product does not mint Member")
    _require((primal * other.dual()).carrier is None, "the swapped mixed product does not take the default carrier")
    both = primal.dual() * other.dual()
    _require(both.carrier == primal.dual().carrier and type(both.member(column)) is Covector, "two dual spaces do not keep the dual carrier")
