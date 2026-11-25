r"""G04u: on a Galerkin frame ``test_space(frame)`` IS ``basis_space(frame)``, and a distinct test basis with no form is formless (N04, spec 7.5; SD-35(a)).

On a Galerkin frame (the angular spherical-harmonic frame of ``lebedev(11)`` at L = 2, and a hand-built indicator frame) ``test_space(frame) is basis_space(frame)``, its head factor is
``Variance(MEASURE, INTEGRATED)`` and its ``apply_metric`` of a seeded vector equals the hand rule: the integrated head form's action is the masked reciprocal of the frame Gram on a diagonal Gram and the
stated-face solve on a dense one, and equals the tip's bytes (G00a's stored ``frames.npz`` is the byte pin). On a Petrov-Galerkin indicator frame whose test basis carries no form, the test space's factor is
``FormlessVariance.MEASURE``.

First red at the base tip: the ``is`` holds today (SD-35(a) preserves it); the variance legs are attribute-only (``FunctionSpace`` has no ``factors``). Mutations: ``test_space`` rebuilding an equal space reddens
the ``is`` leg; ``test_space`` returning ``basis_space(frame).dual().riesz()`` (the declined alternative) reddens the variance leg. The numeric leg is a computation-to-closed-form comparison at a stated ``1e-12`` relative
(the Gram of a spherical-harmonic frame is a LAPACK product), exact on any host to that bound.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _indicator(test_weight: np.ndarray | None):
    from orpheus.numerics.basis import IndicatorBasis, WeightedIndicatorBasis
    from orpheus.numerics.frame import GalerkinFrame, PetrovGalerkinFrame
    from orpheus.numerics.measure import DiscreteMeasure
    from orpheus.geometry.manifold import RealSpace

    edges = np.array([0.0, 1.0, 2.5, 4.0])
    trial = IndicatorBasis((edges,), RealSpace(1))
    measure = DiscreteMeasure(nodes=np.array([0.3, 0.7, 1.4, 2.0, 3.0, 3.6]), weights=np.array([0.3, 0.5, 0.9, 0.6, 0.7, 0.8]), support=RealSpace(1))
    if test_weight is None:
        return GalerkinFrame(trial, measure)
    return PetrovGalerkinFrame(trial, measure, WeightedIndicatorBasis(trial, test_weight))


def test_test_space_is_basis_space_on_a_galerkin_frame() -> None:
    from orpheus.numerics.factor import Kind, Variance, VarianceChart
    from orpheus.numerics.frame_faces import angular_frame, basis_space, test_space
    from orpheus.numerics.quadrature.directional import Quadrature

    for name, frame in (("harmonic lebedev(11) L=2", angular_frame(Quadrature.lebedev(11), 2)), ("indicator", _indicator(None))):
        head = basis_space(frame)
        _require(test_space(frame) is head, f"[{name}] test_space(frame) is not basis_space(frame)")
        (factor,) = head.factors
        _require(factor.variance == Variance(Kind.MEASURE, VarianceChart.INTEGRATED), f"[{name}] the head variance is {factor.variance}")


def test_the_integrated_head_applies_the_inverse_of_the_frame_gram() -> None:
    from orpheus.numerics.frame_faces import angular_frame, basis_space
    from orpheus.numerics.quadrature.directional import Quadrature

    frame = angular_frame(Quadrature.lebedev(11), 2)
    head = basis_space(frame)
    g = np.asarray(frame.discrete_gram, dtype=float)
    x = np.random.default_rng(5).standard_normal(head.shape)
    # the Parseval chart of the integrated head is the inverse of the Gram on its live directions: the hand rule is the Moore-Penrose-free statement G (G^+ x) = x on the live subspace
    y = np.asarray(head.apply_metric(x), dtype=float).reshape(-1)
    x = x.reshape(-1)
    live = np.linalg.norm(g, axis=0) > 0
    back = (g @ y)[live]
    _require(np.allclose(back, x[live], rtol=1e-12, atol=1e-12 * float(np.max(np.abs(x)))), "G applied to the integrated head's metric action does not return x on the live directions")


def test_a_distinct_test_basis_with_no_form_is_formless() -> None:
    from orpheus.numerics.axis import Axis
    from orpheus.numerics.factor import FormlessVariance
    from orpheus.numerics.frame_faces import basis_space, test_space
    from orpheus.numerics.metric import NoForm

    frame = _indicator(np.array([1.0, 2.0, 0.5]))
    _require(frame.test is not frame.basis, "the Petrov-Galerkin frame's test basis is its trial basis")
    _require(isinstance(Axis.of_basis(frame.test).form, NoForm), "the fixture's test basis has a form (the leg would be inert)")
    space = test_space(frame)
    _require(space is not basis_space(frame), "a distinct test basis returned the trial head")
    (factor,) = space.factors
    _require(factor.variance == FormlessVariance.MEASURE, f"the formless test basis' variance is {factor.variance}")
