r"""#707 Step A, the numerics rows: a space mints its members (R5), the flattened seam reads the ends, and the Hilbert
adjoint refuses an operator whose ends lie on opposite sides.

Spec: ``scratch/true_ends/stepA/gate_spec.md`` rows NM-1 to NM-3. Every row is RED today for the reason Step A fixes.
The fixture is a two-axis diagonal-metric space ``V`` (weights 0.5, 2, 4 on one axis, 1 and 3 on the other) and its
dual ``V*``: the smallest spaces on which a chart mistake (a stray or missing metric factor) is not invisible.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest

from orpheus.numerics.axis import Axis
from orpheus.numerics.flat_operator import FlattenedOperator
from orpheus.geometry.manifold import RealSpace
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import LinearOperator, MissingAdjoint
from orpheus.numerics.space import Covector, FunctionSpace
from tests.gates._true_ends import H_REFUSAL_FRAGMENT

pytestmark = pytest.mark.foundation



def _axis(label: str, weights: np.ndarray) -> Axis:
    nodes = np.arange(weights.size, dtype=float)
    return Axis.of_measure(DiscreteMeasure(nodes=nodes, weights=weights, support=RealSpace(1)), label)


def _space() -> FunctionSpace:
    return FunctionSpace.of_axes(_axis("a", np.array([0.5, 2.0, 4.0])), _axis("b", np.array([1.0, 3.0])))


class _CrossSide(LinearOperator):
    """A V → V* double honouring the operator contract: ``apply`` is a fixed dense matrix on coordinates."""

    def __init__(self, V: FunctionSpace) -> None:
        self._V = V
        self._A = np.arange(36.0).reshape(6, 6) / 7.0 + np.eye(6)

    @property
    def domain(self):
        return self._V

    @property
    def codomain(self):
        return self._V.dual()

    @property
    def is_adjointable(self) -> bool:
        return True

    def apply(self, x):
        return (self._A @ np.asarray(x, float).ravel()).reshape(self._V.shape)

    def apply_transpose(self, y):
        return (self._A.T @ np.asarray(y, float).ravel()).reshape(self._V.shape)

    def inverse(self) -> "_CrossSideInverse":
        """One inverse instance per operator (the behavioural NM-2 row spies on it)."""
        if not hasattr(self, "_inv"):
            self._inv = _CrossSideInverse(self)
        return self._inv


class _CrossSideInverse(LinearOperator):
    """``A⁻¹ : V* → V``: reads the coordinates of a member of V* and solves."""

    def __init__(self, forward: _CrossSide) -> None:
        self._f = forward

    @property
    def domain(self):
        return self._f.codomain

    @property
    def codomain(self):
        return self._f.domain

    def apply(self, y):
        coords = self.domain.coordinates_of(y)
        return np.linalg.solve(self._f._A, np.asarray(coords, float).ravel()).reshape(self._f._V.shape)


def test_nm1_a_space_mints_its_members() -> None:
    """NM-1 RED (R5). ``V.member(c)`` is the inverse of V's chart: its coordinates are ``c``, its space is V; and the
    class follows the space: ``V.dual().member(c)`` is a ``Covector`` (the integrated data side, ``space.py:865``).

    First red today: ``FunctionSpace`` has no ``member``.
    """
    V = _space()
    c = np.arange(6.0).reshape(V.shape)
    x = V.member(c)
    assert x.space is V and np.array_equal(x.values, c)
    y = V.dual().member(c)
    assert isinstance(y, Covector) and y.space is V.dual() and np.array_equal(y.values, c)


@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_nm2_the_flattened_seam_has_no_template() -> None:
    """NM-2 RED (design v2 item 6, census A: 19 448 production flux → sweep calls through the template).

    ``FlattenedOperator(inner)`` takes no template: ``apply`` unflattens through ``inner.domain.member``,
    ``apply_transpose`` through the codomain's, ``inverse()`` through the inverse's own ends. First red today: the
    constructor's ``template`` parameter (``flat_operator.py:115``).
    """
    assert "template" not in inspect.signature(FlattenedOperator.__init__).parameters


@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_nm2_the_flattened_inverse_hands_its_inner_a_member_of_the_inverse_domain() -> None:
    """NM-2 RED, behavioural. On a V → V* operator, ``FlattenedOperator(A).inverse()`` hands ``A.inverse()`` a member of
    V*, never of V. First red today: the constructor needs a template, so the template-free form does not exist."""
    V = _space()
    A = _CrossSide(V)
    flat = FlattenedOperator(A)  # type: ignore[call-arg]
    received = []
    inner_inv = A.inverse()
    original = inner_inv.apply

    def spy(x):
        received.append(x)
        return original(x)

    inner_inv.apply = spy  # type: ignore[method-assign]
    flat.inverse().apply(np.ones(V.n_coordinates))
    assert received and received[0].space is V.dual()


@pytest.mark.xfail(strict=True, reason="#707 Step A: red by design until the true-ends carve (ontology chart W1 to W4); the flip-edit removes this marker")
def test_nm3_hilbert_adjoint_refuses_across_sides() -> None:
    """NM-3 RED. ``.H`` on an operator whose ends lie on opposite sides (V → V*) raises ``MissingAdjoint`` naming
    ``DualPairing`` (which owns the adjoint there); ``.dual()`` stays the metric-free arrow V** = V → V*.

    First red today: ``AdjointOperator`` builds both Riesz legs and returns (``operator.py:1672``).
    """
    A = _CrossSide(_space())
    assert A.dual().domain is A.codomain.dual() and A.dual().codomain is A.domain.dual()
    with pytest.raises(MissingAdjoint, match=H_REFUSAL_FRAGMENT):
        A.H


def test_nm3_hilbert_adjoint_stays_on_one_side() -> None:
    """NM-3 GREEN control. A same-side operator V → V keeps ``.H = G⁻¹ Aᵀ G`` (2 of 18 production sites)."""
    V = _space()

    class _Same(_CrossSide):
        @property
        def codomain(self):
            return self._V

    A = _Same(V)
    G = np.asarray(V.metric_diagonal(), float).ravel()
    y = np.linspace(-1.0, 2.0, 6).reshape(V.shape)
    want = (A._A.T @ (G * y.ravel())) / G
    np.testing.assert_allclose(np.asarray(A.H.apply(y)).ravel(), want, rtol=1e-14)
