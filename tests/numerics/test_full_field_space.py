r"""Identity + structure contract for :class:`FullFieldSpace` (Wave O / O.2b R5).

The composite direct-sum space :math:`V_{\rm bulk}\oplus V_{\rm trace}`. These
tests pin the **identity** semantics (``(name, shape)``, block spaces as
``compare=False`` leaf metadata) and the :meth:`from_blocks` shape derivation —
the contract the :class:`~orpheus.numerics.operator.OperatorSum` composition
guard relies on (every SN composite operand reports the SAME composite domain).

The per-block *metric* semantics (``apply_metric`` / ``apply_inverse_metric`` /
``inner_product`` on a real bulk :math:`\oplus` boundary
:class:`~orpheus.transport.timed_full_field.TimedFullField`) are exercised at
the SN layer in
``tests/sn/operators/test_g_adjoint_reciprocity.py`` and against the dense-probe
oracle ``derivations/diagnostics/diag_p42_adjoint_oracle.py`` — kept there to
avoid importing the transport layer into a numerics-layer test.
"""

from __future__ import annotations

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
def test_identity_is_name_and_shape_only():
    # Same total dimension but DIFFERENT block metrics → still equal: the block
    # spaces are compare=False leaf metadata, identity is (name, shape).
    a = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1), weights=True), _trace(8, weights=True))
    b = FullFieldSpace.from_blocks(_bulk((4, 2, 3, 1), weights=False), _trace(8, weights=False))
    assert a == b
    assert hash(a) == hash(b)


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
