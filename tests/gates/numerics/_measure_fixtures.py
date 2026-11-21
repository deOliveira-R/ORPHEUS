"""Test-side fixtures of the numerics measure vocabulary: the values and factories that left production in wave W0, unit ``dead``.

Five definitions moved here from ``orpheus/numerics`` because no production path read them (the dead census, `dead/reviewed.json`
rows 54, 55, 56, 60, 61): the Hermite generating measure :data:`HERMITE`, the generalised Laguerre factory :func:`laguerre`, the
energy-axis manifold :data:`ENERGY`, the midpoint-rule measure :func:`equispaced` and the SymPy-text constructor :func:`from_srepr`.
A moved definition keeps its name and has this one home; a moved method became a free function taking its receiver's place as the
first argument or, for a classmethod, naming its class in the body (decision D2). The helpers a moved body calls and no row rules
(``_laguerre``, ``hermite_recurrence``) stay in production and are imported from it.
"""

from __future__ import annotations

import numpy as np

from orpheus.geometry.manifold import REAL_LINE, Interval
from orpheus.numerics.generating_measure import GeneratingMeasure, _laguerre
from orpheus.numerics.index_set import EnergyGroups
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.mesh_free_function import Symbolic
from orpheus.numerics.recurrence import hermite_recurrence

#: The energy axis of unspecified group count, the support a multigroup measure declares before it knows ``ng``.
ENERGY = EnergyGroups()

#: Weight :math:`w(x) = e^{-x^2}` on :math:`\mathbb{R}`.
HERMITE = GeneratingMeasure(
    name="hermite",
    support=REAL_LINE,
    exact_recurrence=hermite_recurrence,
)


def laguerre(a: float = 0.0) -> GeneratingMeasure:
    """The generalised Laguerre measure for exponent ``a``: the production ``_laguerre`` at ``float(a)``.

    Cached per ``float(a)`` by ``_laguerre``, so a repeated call returns the same measure.
    """
    return _laguerre(float(a))


def equispaced(a: float, b: float, n: int) -> DiscreteMeasure:
    r"""Equispaced (midpoint-rule) measure on :math:`[a, b]`: ``n`` nodes at the midpoints of an equal partition, each of weight
    :math:`(b - a)/n`, on ``support=Interval(a, b)`` with ``degree_of_exactness=1`` against Lebesgue measure on that interval.

    The same nodes read as a rule on the circle are the periodic trapezoid, a different claim about a different reference, which
    belongs to a circle rule; this one is the interval's midpoint rule. Refused: ``n < 1`` and ``a >= b``.
    """
    if n < 1:
        raise ValueError(f"equispaced requires n >= 1, got n={n}")
    if not (a < b):
        raise ValueError(f"equispaced requires a < b, got a={a}, b={b}")
    from orpheus.numerics.exactness import ExactnessClaim, OrthogonalSystem, UniformMeasure

    h = (b - a) / n
    nodes = a + (np.arange(n) + 0.5) * h
    weights = np.full(n, h)
    support = Interval(a, b)
    return DiscreteMeasure(
        nodes=nodes,
        weights=weights,
        support=support,
        exactness=ExactnessClaim(
            reference=UniformMeasure(support=support, orthogonal_system=OrthogonalSystem.ALGEBRAIC),
            degree=1,
        ),
    )


def from_srepr(srepr: tuple[str, ...], sympy_version: str | None = None) -> Symbolic:
    """Stored text, with the version of the SymPy that wrote it (this process's when omitted).

    Was ``Symbolic.from_srepr``, a classmethod whose class parameter every call site filled with ``Symbolic``: the class is named in
    the body and the parameter is dropped.
    """
    import sympy

    return Symbolic(srepr=tuple(srepr), sympy_version=sympy.__version__ if sympy_version is None else sympy_version)
