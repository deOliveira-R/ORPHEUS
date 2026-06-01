"""Foundation tests for :mod:`orpheus.numerics.space`.

Tests the FunctionSpace primitive shipped in Phase B / Issue 9.6 of
the SN reshape campaign — equality / hash semantics on
``(name, shape)``, repr formatting, weighted vs Euclidean inner
product, norm consistency, and the pre-populated common-space
factories.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.space import (
    FunctionSpace,
    angular_flux_space,
    scalar_flux_space,
)


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
def test_function_space_equality_independent_of_weights():
    """Different inner-product weights does NOT break identity."""
    w = np.arange(4.0)
    a = FunctionSpace(name="psi", shape=(10, 4))
    b = FunctionSpace(name="psi", shape=(10, 4), inner_product_weights=w)
    assert a == b
    # And they hash the same.
    assert hash(a) == hash(b)


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


# ─────────────────────────────────────────────────────────────────────
# Pre-populated common-space factories
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_angular_flux_space_shape():
    sp = angular_flux_space(n_cells=10, n_ordinates=8, n_groups=2)
    assert sp.shape == (10, 8, 2)
    assert sp.name == "angular_flux"
    assert sp.inner_product_weights is None


@pytest.mark.foundation
def test_angular_flux_space_with_quadrature_weights():
    n_ord = 8
    qw = np.linspace(0.1, 0.9, n_ord)
    sp = angular_flux_space(
        n_cells=10, n_ordinates=n_ord, n_groups=2,
        quadrature_weights=qw,
    )
    assert sp.inner_product_weights is not None
    assert sp.inner_product_weights.shape == (1, n_ord, 1)


@pytest.mark.foundation
def test_angular_flux_space_quadrature_weights_shape_check():
    with pytest.raises(ValueError, match="quadrature_weights"):
        angular_flux_space(
            n_cells=10, n_ordinates=8, n_groups=2,
            quadrature_weights=np.ones(7),  # wrong length
        )


@pytest.mark.foundation
def test_scalar_flux_space_shape():
    sp = scalar_flux_space(n_cells=10, n_groups=2)
    assert sp.shape == (10, 2)
    assert sp.name == "scalar_flux"
    assert sp.inner_product_weights is None


# NOTE: the directional ``boundary_trace_space()`` factory was retired
# in the #205/#201 trace-space unification — the whole-boundary trace is
# now the concrete :class:`~orpheus.numerics.spaces.trace_space.TraceSpace`
# (see ``tests/numerics/test_trace_space.py``), with inflow/outflow as
# selectors rather than a directional space tag.
