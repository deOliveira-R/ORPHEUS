"""Constructors the name-built fixtures of the numerics gates share (N04, spec 6.5 "Tests, the file set"): not a test module.

A fixture whose form holds a zero weight, or a dense form with a declared kernel, is a space over a singular form, which a ``Factor`` admits
only on the quotient by the form's kernel: the whole index set ``S``, modulo ``K`` the form's own kernel (the zero indices of a diagonal form,
the declared kernel of a dense one). ``declared_mod_kernel`` is that one spelling, so each such fixture names its subquotient in its call.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from orpheus.numerics.axis import Axis, BasisKind, Declared
from orpheus.numerics.factor import Factor, Kind, Subquotient, Variance, VarianceChart
from orpheus.numerics.metric import DiagonalMetric, NoForm
from orpheus.numerics.factor import Kind as _Kind
from orpheus.numerics.space import DualSpace, TensorSpace
from orpheus.numerics.declaration import Coord


def declared_mod_kernel(name: str, shape: tuple[int, ...], variance: Any, *, form: Any) -> TensorSpace:
    """``TensorSpace.declared(name, shape, variance, form=form)`` on the subquotient ``Coord(n, all) / form.kernel``."""
    axis = Axis(name, tuple(int(d) for d in shape), origin=Declared(BasisKind.NODAL, form))
    n = int(np.prod(axis.shape))
    return TensorSpace.of_axes(Factor(axis, variance, Subquotient(Coord(n, tuple(range(n))), form.kernel)))


def factor_mod_kernel(axis: Axis, variance: Any = None) -> Factor:
    """The factor of ``axis`` on the subquotient ``Coord(n, all) / axis.form.kernel`` (``(MEASURE, DENSITY)`` unless ``variance`` is given)."""
    n = int(np.prod(axis.shape))
    v = Variance(Kind.MEASURE, VarianceChart.DENSITY) if variance is None else variance
    return Factor(axis, v, Subquotient(Coord(n, tuple(range(n))), axis.form.kernel))


def whole_on_mod_kernel(axis: Axis, kind: Kind, chart: VarianceChart) -> Factor:
    """``Factor.whole_on(axis, kind, chart)`` on the quotient by the axis form's kernel where the axis has one (a form with a zero weight)."""
    if isinstance(axis.form, NoForm):
        return Factor.whole_on(axis, kind, chart)
    n = int(np.prod(axis.shape))
    return Factor(axis, Variance(kind, chart), Subquotient(Coord(n, tuple(range(n))), axis.form.kernel))


def leading_aligned(weights: Any, shape: tuple[int, ...]) -> DiagonalMetric:
    """The diagonal form of ``weights`` read leading-aligned over ``shape``: a proper prefix of the shape (or any broadcastable array) broadcast to the whole shape.

    The weights of the retired ``inner_product_weights=`` slot were broadcast against the leading axes of an element; a factor's form lives over exactly
    the axis shape, so the broadcast is stated here."""
    w = np.asarray(weights, dtype=float)
    return DiagonalMetric(np.broadcast_to(w.reshape(w.shape + (1,) * (len(shape) - w.ndim)), tuple(int(d) for d in shape)))


def legacy_leaf(name: str, shape: tuple[int, ...]) -> DualSpace:
    """A legacy leaf space (no named factors) over ``shape``: the one arm of ``FunctionSpace`` a test can build that has no factors to look up (a ``DualSpace``)."""
    primal = TensorSpace.declared(name, shape, Variance(_Kind.FUNCTION, VarianceChart.DENSITY))
    return DualSpace(name=name, shape=tuple(int(d) for d in shape), primal=primal)
