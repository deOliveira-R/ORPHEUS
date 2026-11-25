r"""Constructors the N04 gates share (spec 7.5). Not a test module: every function RETURNS data and none asserts.

``python -O`` strips a bare ``assert`` of a helper, so the collected test files assert. Every name that exists only after the unit
(``TensorSpace``, ``CarrierSpec``) is read from its module AT CALL TIME, inside each function, so that on a tree without it a gate row fails in
its body (a behavioural red a reader can interpret) and not at collection.
"""

from __future__ import annotations

from typing import Any

import numpy as np

# The angular, energy and spatial weights of the small slab-like triple every N04 gate builds (non-uniform, so that two legs applied in two
# orders differ in the last bit and a Kronecker walk in the wrong order is visible).
ANGULAR_W = (0.7, 1.3, 0.45, 2.1)
SPATIAL_W = (1.1, 0.35, 2.7)
SIGNED_A = (1.3, -0.7, 2.1)
SIGNED_B = (0.3, -1.1, 0.9, 0.0)


def tensor_space() -> Any:
    from orpheus.numerics.space import TensorSpace

    return TensorSpace


def diag_axis(weights: tuple[float, ...] | np.ndarray, label: str) -> Any:
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.metric import DiagonalMetric

    w = np.asarray(weights, dtype=float)
    return Axis(label, (len(w),), origin=Declared(BasisKind.NODAL, DiagonalMetric(w)))


def formless_axis(n: int, label: str) -> Any:
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.metric import NoForm

    return Axis(label, (n,), origin=Declared(BasisKind.NODAL, NoForm()))


def density_factor(axis: Any, *, sk: Any = None) -> Any:
    """The ``(MEASURE, DENSITY)`` factor of ``axis`` (whole unless ``sk`` is given)."""
    from orpheus.numerics import factor as F

    n = int(axis.shape[0])
    return F.Factor(axis, F.Variance(F.Kind.MEASURE, F.VarianceChart.DENSITY), F.Subquotient.whole(n) if sk is None else sk)


def formless_factor(axis: Any, variance: Any = None) -> Any:
    """The factor of a formless axis; ``FormlessVariance.MEASURE`` unless ``variance`` is given."""
    from orpheus.numerics import factor as F

    return F.Factor(axis, F.FormlessVariance.MEASURE if variance is None else variance, F.Subquotient.whole(int(axis.shape[0])))


def triple_factors(angular_w: tuple[float, ...] = ANGULAR_W, ng: int = 2, spatial_w: tuple[float, ...] = SPATIAL_W) -> tuple[Any, Any, Any]:
    """``(angular, energy, spatial)`` factors in the slab's slot order: angular and spatial ``(MEASURE, DENSITY)``, energy ``FormlessVariance.MEASURE``."""
    return (
        density_factor(diag_axis(angular_w, "angular")),
        formless_factor(formless_axis(ng, "energy")),
        density_factor(diag_axis(spatial_w, "spatial")),
    )


def triple_space(angular_w: tuple[float, ...] = ANGULAR_W) -> Any:
    return tensor_space().of_axes(*triple_factors(angular_w))


def signed_pair_space() -> Any:
    """Two formed slots whose forms are signed: ``Diag(SIGNED_A)`` and ``Diag(SIGNED_B)`` (a zero weight, so the second takes the quotient by its kernel)."""
    from orpheus.numerics import factor as F
    from orpheus.numerics.declaration import Coord

    a = density_factor(diag_axis(SIGNED_A, "a"))
    kernel_free = Coord(4, (0, 1, 2, 3))
    b = density_factor(diag_axis(SIGNED_B, "b"), sk=F.Subquotient(kernel_free, Coord(4, (3,))))
    return tensor_space().of_axes(a, b)


def slot_order_flat(space: Any, x: np.ndarray, labels: tuple[str, ...], *, reverse: bool = False) -> np.ndarray:
    """A HAND spelling of ``flat_on`` for diagonal forms: one broadcast product per formed slot in slot order (or reversed), written from the definition.

    It reads each factor's form weights and applies ``y * w`` along the slot's array axis; a ``NoForm`` slot is skipped.
    """
    from orpheus.numerics.metric import DiagonalMetric

    slots = [f for f in space.factors if f.axis.label in labels]
    if reverse:
        slots = slots[::-1]
    y = np.array(x, dtype=float)
    for f in slots:
        form = f.form
        if not isinstance(form, DiagonalMetric):
            continue
        position = [g.axis.label for g in space.factors].index(f.axis.label)
        shape = [1] * y.ndim
        shape[position] = -1
        y = y * np.asarray(form.weights, dtype=float).reshape(shape)
    return y
