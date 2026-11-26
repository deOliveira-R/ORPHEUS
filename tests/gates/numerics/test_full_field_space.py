r"""Identity + structure contract for :class:`DirectSumSpace` (Wave O / O.2b R5).

The composite direct-sum space :math:`V_{\rm bulk}\oplus V_{\rm trace}`. These
tests pin the **identity** semantics (the members' keys: a composite is its
blocks, and an axes-less block's identity carries its metric) and the
:meth:`of_pairs` shape derivation —
the contract the :class:`~orpheus.numerics.operator.OperatorSum` composition
guard relies on (every SN composite operand reports the SAME composite domain).

The per-block *metric* semantics (``apply_metric`` / ``apply_inverse_metric`` /
``inner_product`` on a real bulk :math:`\oplus` boundary
:class:`~orpheus.transport.timed_full_field.TimedFullField`) are exercised at
the SN layer in
``tests/gates/sn/operators/test_g_adjoint_reciprocity.py`` — kept there to avoid
importing the transport layer into a numerics-layer test. (A dense-probe
oracle, ``diag_p42_adjoint_oracle.py``, checked them at landing; it was not
preserved in git.)
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from orpheus.numerics.space import FunctionSpace
from orpheus.numerics.direct_sum import DirectSumSpace
from orpheus.numerics.factor import Kind, Variance, VarianceChart
from orpheus.numerics.space import TensorSpace
from orpheus.numerics.metric import DiagonalMetric
from orpheus.numerics.direct_sum import BULK, DirectSumSpace, TRACE


def _bulk(shape=(4, 2, 3, 1), weights=True) -> FunctionSpace:
    w = np.arange(1.0, float(np.prod(shape)) + 1.0).reshape(shape) if weights else None
    return TensorSpace.declared("sn_bulk", shape, Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=None if w is None else DiagonalMetric(w))


def _trace(n=8, weights=True) -> FunctionSpace:
    w = np.linspace(0.1, 1.0, n) if weights else None
    return TensorSpace.declared("angular_trace", (n,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=None if w is None else DiagonalMetric(w))


@pytest.mark.foundation
def test_from_blocks_shape_is_flat_direct_sum():
    space = DirectSumSpace.of_pairs(((BULK, _bulk((4, 2, 3, 1))), (TRACE, _trace(8))))
    # Family prefix + a member-content digest (CS4b S4 — the of_axes rule
    # for direct sums); never pin the full literal (R4).
    assert space.name.startswith("dsum(bulk ⊕ trace)#")
    assert space.shape == (4 * 2 * 3 * 1 + 8,)
    assert space.summand("bulk").shape == (4, 2, 3, 1)
    assert space.summand("trace").shape == (8,)


@pytest.mark.foundation
def test_identity_is_the_members_content():
    # A composite's identity is its blocks' identities, and an axes-less block's identity carries its metric
    # (U1, 2026-10-30): equal blocks give equal composites; blocks that differ only in their metric give two
    # spaces. Until U1 this pair compared EQUAL, the metric-blind (name, shape) identity.
    a = DirectSumSpace.of_pairs(((BULK, _bulk((4, 2, 3, 1), weights=True)), (TRACE, _trace(8, weights=True))))
    same = DirectSumSpace.of_pairs(((BULK, _bulk((4, 2, 3, 1), weights=True)), (TRACE, _trace(8, weights=True))))
    b = DirectSumSpace.of_pairs(((BULK, _bulk((4, 2, 3, 1), weights=False)), (TRACE, _trace(8, weights=False))))
    assert a == same and hash(a) == hash(same)
    assert a != b


@pytest.mark.foundation
def test_inequality_on_different_total_dimension():
    a = DirectSumSpace.of_pairs(((BULK, _bulk((4, 2, 3, 1))), (TRACE, _trace(8))))
    b = DirectSumSpace.of_pairs(((BULK, _bulk((4, 2, 3, 1))), (TRACE, _trace(9))))
    assert a != b


@pytest.mark.foundation
def test_usable_as_dict_key():
    a = DirectSumSpace.of_pairs(((BULK, _bulk()), (TRACE, _trace())))
    b = DirectSumSpace.of_pairs(((BULK, _bulk()), (TRACE, _trace())))  # equal identity
    d = {a: "composite"}
    assert d[b] == "composite"


@pytest.mark.foundation
def test_repr_names_the_space_and_shape():
    space = DirectSumSpace.of_pairs(((BULK, _bulk((4, 2, 3, 1))), (TRACE, _trace(8))))
    r = repr(space)
    assert "DirectSumSpace" in r
    assert "dsum(bulk ⊕ trace)" in r
    assert str(space.shape) in r


@pytest.mark.foundation
def test_composite_own_inner_product_weights_is_refused():
    # The composite carries no own diagonal metric — the metric lives per block — and, as on a tensor space, reading the stored field is refused.
    space = DirectSumSpace.of_pairs(((BULK, _bulk()), (TRACE, _trace())))
    with pytest.raises(AttributeError, match="stores no metric source"):
        space.inner_product_weights


@pytest.mark.foundation
def test_a_composite_block_carrying_a_dense_metric_delegates_it():
    r"""C5 (P7) — verify-don't-assume on the composite delegation: a
    DENSE-metric interior block must flow through the direct-sum pairing.

    Hand-derived exact-binary literals (the dense-metric battery's G):
    interior ⟨x,y⟩_G = 93/4 = 23.25, trace ⟨x,y⟩_w = 2 + 1/2 = 2.5,
    sum compared with ``==``. The three composite verbs share one
    delegation line each (``<block>_space.<verb>(x.<block>.values)``),
    so the pairing face witnesses the delegation for all three.
    """
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.metric import DenseMetric, Face

    g = np.array([[2.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 4.0]])
    interior = TensorSpace.declared("bulk3", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DenseMetric(g, kernel=Coord(3, ()), face=Face.METRIC))
    trace = TensorSpace.declared("tr2", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(np.array([2.0, 0.5])))
    space = DirectSumSpace.of_pairs(((BULK, interior), (TRACE, trace)))
    x = _Composite(_Leaf(np.array([1.0, 2.0, 3.0])), _Leaf(np.array([1.0, 2.0])))
    y = _Composite(_Leaf(np.array([0.5, -1.0, 2.0])), _Leaf(np.array([1.0, 0.5])))
    assert space.inner_product(x, y) == 23.25 + 2.5


@dataclasses.dataclass(frozen=True)
class _Leaf:
    values: np.ndarray


@dataclasses.dataclass(frozen=True)
class _Composite:
    """A minimal carrier double: two frozen leaves and the ``blocks`` and ``with_blocks`` accesses the sum's legs read and rebuild through."""

    interior: _Leaf
    boundary: _Leaf

    @property
    def blocks(self):
        return (self.interior, self.boundary)

    def with_blocks(self, blocks):
        return _Composite(*blocks)


@pytest.mark.parametrize("leg", ["apply_metric", "apply_inverse_metric"])
def test_the_flat_chart_and_the_carrier_agree_bit_for_bit(leg):
    r"""An element of the direct sum has two charts, the typed carrier and the flat coordinates
    ``[interior.ravel(), trace]``, and each metric leg acts on both (#666's flat half, U0). The flat chart is what the
    operator algebra hands the legs (an operator's ``.H``, the kernel gauge's ``gather.H``). The two results agree
    bit for bit, and a trailing payload rides through as ``G (x) I_P``.

    First red: before U0 the flat chart raised ``AttributeError: 'numpy.ndarray' object has no attribute
    'interior'``. Arm: the blocks concatenated trace-first reddens the agreement."""
    space = DirectSumSpace.of_pairs(((BULK, _bulk()), (TRACE, _trace())))
    rng = np.random.default_rng(666)
    bulk, trace = rng.standard_normal((4, 2, 3, 1)), rng.standard_normal(8)
    typed = getattr(space, leg)(_Composite(_Leaf(bulk), _Leaf(trace)))
    flat = getattr(space, leg)(np.concatenate([bulk.ravel(), trace]))
    np.testing.assert_array_equal(flat, np.concatenate([typed.interior.values.ravel(), typed.boundary.values]))
    batch = rng.standard_normal((space.shape[0], 3))
    columns = np.column_stack([getattr(space, leg)(batch[:, j]) for j in range(3)])
    np.testing.assert_array_equal(getattr(space, leg)(batch), columns)


def test_the_flat_chart_is_parsed_at_the_composite():
    """A flat operand of the wrong length is not a coordinate vector of the sum: refused, not sliced short."""
    space = DirectSumSpace.of_pairs(((BULK, _bulk()), (TRACE, _trace())))
    with pytest.raises(ValueError, match="shape"):
        space.apply_metric(np.ones(space.shape[0] - 1))
