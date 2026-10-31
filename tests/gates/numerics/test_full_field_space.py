r"""Identity + structure contract for :class:`FullFieldSpace` (Wave O / O.2b R5).

The composite direct-sum space :math:`V_{\rm bulk}\oplus V_{\rm trace}`. These
tests pin the **identity** semantics (the members' keys: a composite is its
blocks, and an axes-less block's identity carries its metric) and the
:meth:`from_blocks` shape derivation —
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
from orpheus.numerics.spaces.full_field_space import FullFieldSpace


def _bulk(shape=(4, 2, 3, 1), weights=True) -> FunctionSpace:
    w = np.arange(1.0, float(np.prod(shape)) + 1.0).reshape(shape) if weights else None
    return FunctionSpace(name="sn_bulk", shape=shape, inner_product_weights=w)


def _trace(n=8, weights=True) -> FunctionSpace:
    w = np.linspace(0.1, 1.0, n) if weights else None
    return FunctionSpace(name="angular_trace", shape=(n,), inner_product_weights=w)


@pytest.mark.foundation
def test_from_blocks_shape_is_flat_direct_sum():
    space = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1)), _trace(8))
    # Family prefix + a member-content digest (CS4b S4 — the of_axes rule
    # for direct sums); never pin the full literal (R4).
    assert space.name.startswith("full_field#")
    assert space.shape == (4 * 2 * 3 * 1 + 8,)
    assert space.interior_space.shape == (4, 2, 3, 1)
    assert space.trace_space.shape == (8,)


@pytest.mark.foundation
def test_identity_is_the_members_content():
    # A composite's identity is its blocks' identities, and an axes-less block's identity carries its metric
    # (U1, 2026-10-30): equal blocks give equal composites; blocks that differ only in their metric give two
    # spaces. Until U1 this pair compared EQUAL, the metric-blind (name, shape) identity.
    a = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1), weights=True), _trace(8, weights=True))
    same = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1), weights=True), _trace(8, weights=True))
    b = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1), weights=False), _trace(8, weights=False))
    assert a == same and hash(a) == hash(same)
    assert a != b


@pytest.mark.foundation
def test_inequality_on_different_total_dimension():
    a = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1)), _trace(8))
    b = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1)), _trace(9))
    assert a != b


@pytest.mark.foundation
def test_usable_as_dict_key():
    a = FullFieldSpace.from_blocks(_bulk(), _trace())
    b = FullFieldSpace.from_blocks(_bulk(), _trace())  # equal identity
    d = {a: "composite"}
    assert d[b] == "composite"


@pytest.mark.foundation
def test_repr_names_the_space_and_shape():
    space = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1)), _trace(8))
    r = repr(space)
    assert "FullFieldSpace" in r
    assert "full_field" in r
    assert str(space.shape) in r


@pytest.mark.foundation
def test_composite_own_inner_product_weights_stays_none():
    # The composite carries no own diagonal metric — the metric lives per block.
    space = FullFieldSpace.from_blocks(_bulk(), _trace())
    assert space.inner_product_weights is None


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
    from types import SimpleNamespace

    from orpheus.numerics.metric import DenseMetric

    g = np.array([[2.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 4.0]])
    interior = FunctionSpace("bulk3", (3,), metric=DenseMetric(g))
    trace = FunctionSpace(
        "tr2", (2,), inner_product_weights=np.array([2.0, 0.5])
    )
    space = FullFieldSpace.from_blocks(interior, trace)
    x = SimpleNamespace(
        interior=SimpleNamespace(values=np.array([1.0, 2.0, 3.0])),
        boundary=SimpleNamespace(values=np.array([1.0, 2.0])),
    )
    y = SimpleNamespace(
        interior=SimpleNamespace(values=np.array([0.5, -1.0, 2.0])),
        boundary=SimpleNamespace(values=np.array([1.0, 0.5])),
    )
    assert space.inner_product(x, y) == 23.25 + 2.5


@dataclasses.dataclass(frozen=True)
class _Leaf:
    values: np.ndarray


@dataclasses.dataclass(frozen=True)
class _Composite:
    """A minimal carrier double: two frozen leaves and the ``_recombine`` hook the composite's legs rebuild through."""

    interior: _Leaf
    boundary: _Leaf

    def _recombine(self, *, interior, boundary):
        return _Composite(interior, boundary)


@pytest.mark.parametrize("leg", ["apply_metric", "apply_inverse_metric"])
def test_the_flat_chart_and_the_carrier_agree_bit_for_bit(leg):
    r"""An element of the direct sum has two charts, the typed carrier and the flat coordinates
    ``[interior.ravel(), trace]``, and each metric leg acts on both (#666's flat half, U0). The flat chart is what the
    operator algebra hands the legs (an operator's ``.H``, the kernel gauge's ``gather.H``). The two results agree
    bit for bit, and a trailing payload rides through as ``G (x) I_P``.

    First red: before U0 the flat chart raised ``AttributeError: 'numpy.ndarray' object has no attribute
    'interior'``. Arm: the blocks concatenated trace-first reddens the agreement."""
    space = FullFieldSpace.from_blocks(_bulk(), _trace())
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
    space = FullFieldSpace.from_blocks(_bulk(), _trace())
    with pytest.raises(ValueError, match="shape"):
        space.apply_metric(np.ones(space.shape[0] - 1))
