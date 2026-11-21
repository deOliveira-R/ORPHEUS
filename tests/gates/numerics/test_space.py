"""Foundation tests for :mod:`orpheus.numerics.space`.

Tests the FunctionSpace primitive shipped in Phase B / Issue 9.6 of
the SN reshape campaign — equality / hash semantics on
``(name, shape)``, repr formatting, weighted vs Euclidean inner
product, norm consistency.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.space import FunctionSpace


# ─────────────────────────────────────────────────────────────────────
# Equality / hash semantics
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_function_space_equality_on_name_and_shape():
    """Two FunctionSpaces with same name and shape compare equal."""
    a = FunctionSpace(name="phi", shape=(10, 4))
    b = FunctionSpace(name="phi", shape=(10, 4))
    assert a == b


@pytest.mark.foundation
def test_function_space_equality_carries_the_weights():
    """The inner-product weights are part of a space's identity (the
    space-membership unit, U1): two spaces of one name and shape with
    different weights are two spaces, and equal weights give equal,
    hash-equal spaces.

    The weights are leading-aligned on the ``(10, 4)`` element (one weight per
    row). Until U1b the fixture carried ``(4,)`` weights, which no element of
    the space could pair with (the leading-aligned ``(4, 1)`` does not
    broadcast to ``(10, 4)``); the space now checks its weights when built."""
    w = np.arange(1.0, 11.0)
    a = FunctionSpace(name="psi", shape=(10, 4))
    b = FunctionSpace(name="psi", shape=(10, 4), inner_product_weights=w)
    assert a != b
    twin = FunctionSpace(name="psi", shape=(10, 4), inner_product_weights=w.copy())
    assert b == twin
    assert hash(b) == hash(twin)


@pytest.mark.foundation
def test_function_space_inequality_on_different_name():
    a = FunctionSpace(name="phi", shape=(10,))
    b = FunctionSpace(name="psi", shape=(10,))
    assert a != b


@pytest.mark.foundation
def test_function_space_inequality_on_different_shape():
    a = FunctionSpace(name="phi", shape=(10,))
    b = FunctionSpace(name="phi", shape=(20,))
    assert a != b


@pytest.mark.foundation
def test_function_space_hash_stable():
    """Hash is consistent with equality."""
    a = FunctionSpace(name="phi", shape=(10, 4))
    b = FunctionSpace(name="phi", shape=(10, 4))
    assert hash(a) == hash(b)


@pytest.mark.foundation
def test_function_space_usable_as_dict_key():
    """FunctionSpace must be hashable for use as a domain/range tag."""
    sp = FunctionSpace(name="phi", shape=(10,))
    d = {sp: "value"}
    assert d[FunctionSpace(name="phi", shape=(10,))] == "value"


# ─────────────────────────────────────────────────────────────────────
# Repr
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_function_space_repr_format():
    sp = FunctionSpace(name="angular_flux", shape=(50, 8, 2))
    r = repr(sp)
    assert "FunctionSpace" in r
    assert "'angular_flux'" in r
    assert "(50, 8, 2)" in r


# ─────────────────────────────────────────────────────────────────────
# Inner product / norm
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_inner_product_euclidean_default():
    """Without weights, inner_product is the Euclidean sum."""
    sp = FunctionSpace(name="phi", shape=(4,))
    x = np.array([1.0, 2.0, 3.0, 4.0])
    y = np.array([1.0, 1.0, 1.0, 1.0])
    assert sp.inner_product(x, y) == pytest.approx(10.0)


@pytest.mark.foundation
def test_inner_product_with_weights():
    """With diagonal weights, it's the weighted sum w_i x_i y_i."""
    w = np.array([1.0, 2.0, 3.0, 4.0])
    sp = FunctionSpace(
        name="phi", shape=(4,), inner_product_weights=w,
    )
    x = np.array([1.0, 1.0, 1.0, 1.0])
    y = np.array([1.0, 1.0, 1.0, 1.0])
    # Sum w_i = 1 + 2 + 3 + 4 = 10.
    assert sp.inner_product(x, y) == pytest.approx(10.0)


@pytest.mark.foundation
def test_inner_product_broadcast_over_extra_axes():
    """Weights along the ordinate axis broadcast over extra axes."""
    n_ord = 4
    n_groups = 2
    n_cells = 3
    w_ord = np.array([1.0, 2.0, 3.0, 4.0])
    w_broadcast = w_ord.reshape(1, n_ord, 1)
    sp = FunctionSpace(
        name="psi",
        shape=(n_cells, n_ord, n_groups),
        inner_product_weights=w_broadcast,
    )
    x = np.ones((n_cells, n_ord, n_groups))
    y = np.ones((n_cells, n_ord, n_groups))
    # Each (cell, group) pair contributes Sum w_i = 10; n_cells * n_groups
    # = 6 such pairs.
    assert sp.inner_product(x, y) == pytest.approx(60.0)


@pytest.mark.foundation
def test_norm_consistent_with_inner_product():
    """norm(x) == sqrt(inner_product(x, x))."""
    sp = FunctionSpace(name="phi", shape=(5,))
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    expected = np.sqrt(1 + 4 + 9 + 16 + 25)
    assert sp.norm(x) == pytest.approx(expected)


# NOTE: the directional ``boundary_trace_space()`` factory was retired
# in the #205/#201 trace-space unification — the whole-boundary trace is
# now the concrete :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
# (see ``tests/gates/numerics/test_trace_space.py``), with inflow/outflow as
# selectors rather than a directional space tag.
