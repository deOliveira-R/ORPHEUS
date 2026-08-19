r"""Foundation tests for the L1 :class:`~orpheus.numerics.field.Field` ABC.

Pins the algebraic contract that every Depth-B typed field inherits
unchanged:

* Same-class + same-space addition / subtraction / negation return
  the same concrete subclass.
* Cross-class arithmetic raises ``TypeError`` (the Layer 1 gate —
  the runtime dimensional defense; runs in both ``pytest`` and
  ``python -O -m pytest``). Under View-G (issues #205 / #207) units
  live on the field role-leaf, not the space, so class identity *is*
  units identity.
* Cross-space arithmetic raises ``ValueError`` even within the same
  class (different shapes or names).
* Scalar multiplication / division / unary negation preserve class
  identity and space.
* :meth:`~Field.inner_product` reads the space's metric (Euclidean
  default; weighted when ``inner_product_weights`` is set).
* Direct instantiation of the abstract :class:`Field` is rejected.
* ``values.shape`` is validated against ``space.shape`` at
  construction.

These tests use a test-local concrete subclass ``_DummyField`` because
:class:`Field` is abstract by design. The concrete production fields
(``ScalarFlux``, ``AngularFlux``, …) are migrated to inherit from
:class:`Field` in later Depth-B steps (D-D through D-H) and carry
their own tests in ``tests/transport/fields/``.

See ``.claude/plans/depth_b_field_on_function_space.md`` §3.2 for the
ABC spec and §7.1 for the test inventory; the View-G dimensional-
enforcement story (units live on the field role-leaf, not the space)
is in the field-vocabulary plan (issues #205 / #201).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from orpheus.numerics.field import Field
from orpheus.numerics.space import FunctionSpace


# ─────────────────────────────────────────────────────────────────────
# Test fixtures: concrete subclasses.
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, eq=False, kw_only=True)
class _DummyField(Field):
    """Minimal concrete subclass — no fields beyond Field's (values, space).

    Used to exercise the ABC's algebra in isolation, without dragging
    in any L2 transport field's added machinery.
    """


@dataclass(frozen=True, eq=False, kw_only=True)
class _OtherDummyField(Field):
    """A distinct concrete subclass with the same field signature as
    ``_DummyField``.

    The class identity invariant rejects cross-class arithmetic even
    when shape and units agree.
    """


@dataclass(frozen=True, eq=False, kw_only=True)
class _RichField(Field):
    """A concrete subclass that adds a non-defaulted extra field.

    Verifies that ``dataclasses.replace`` correctly propagates extra
    fields through the dunder algebra without manual override. The
    extra ``tag: str`` simulates the mesh / boundary / history fields
    that the L2 types will carry.
    """

    tag: str


# ─────────────────────────────────────────────────────────────────────
# Construction-time validation
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_field_is_abstract_direct_instantiation_rejected():
    """The ABC itself MUST refuse direct instantiation."""
    sp = FunctionSpace(name="dummy", shape=(3,))
    with pytest.raises(TypeError, match="abstract"):
        Field(values=np.zeros(3), space=sp)


@pytest.mark.foundation
def test_field_post_init_validates_shape_mismatch():
    sp = FunctionSpace(name="dummy", shape=(3, 2))
    with pytest.raises(ValueError, match="does not match space.shape"):
        _DummyField(values=np.zeros(5), space=sp)


@pytest.mark.foundation
def test_field_post_init_accepts_matching_shape():
    sp = FunctionSpace(name="dummy", shape=(2, 3))
    f = _DummyField(values=np.zeros((2, 3)), space=sp)
    assert f.values.shape == (2, 3)
    assert f.space is sp


# ─────────────────────────────────────────────────────────────────────
# Same-class, same-space algebra
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_addition_returns_same_class_and_space():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.array([1.0, 2.0, 3.0]), space=sp)
    b = _DummyField(values=np.array([10.0, 20.0, 30.0]), space=sp)
    c = a + b
    assert type(c) is _DummyField
    assert c.space is sp
    np.testing.assert_array_equal(c.values, [11.0, 22.0, 33.0])


@pytest.mark.foundation
def test_subtraction_returns_same_class_and_space():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.array([10.0, 20.0, 30.0]), space=sp)
    b = _DummyField(values=np.array([1.0, 2.0, 3.0]), space=sp)
    c = a - b
    assert type(c) is _DummyField
    np.testing.assert_array_equal(c.values, [9.0, 18.0, 27.0])


@pytest.mark.foundation
def test_negation_returns_same_class():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.array([1.0, -2.0, 3.0]), space=sp)
    c = -a
    assert type(c) is _DummyField
    np.testing.assert_array_equal(c.values, [-1.0, 2.0, -3.0])


@pytest.mark.foundation
def test_scalar_multiplication_left_and_right():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.array([1.0, 2.0, 3.0]), space=sp)
    left = 2.5 * a
    right = a * 2.5
    assert type(left) is _DummyField
    assert type(right) is _DummyField
    np.testing.assert_array_equal(left.values, [2.5, 5.0, 7.5])
    np.testing.assert_array_equal(right.values, [2.5, 5.0, 7.5])


@pytest.mark.foundation
def test_scalar_division():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.array([2.0, 4.0, 8.0]), space=sp)
    c = a / 2.0
    assert type(c) is _DummyField
    np.testing.assert_array_equal(c.values, [1.0, 2.0, 4.0])


# ─────────────────────────────────────────────────────────────────────
# Cross-class rejection (Layer 1 gate)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_cross_class_addition_rejected():
    """Layer 1: ``type(self) is not type(other)`` raises TypeError.

    Even when ``space``, ``values.shape``, and ``units`` agree, two
    distinct Field subclasses MUST NOT add. Per plan §3.2:
    "Cross-class arithmetic (even with matching units) requires an
    explicit named composition."
    """
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.zeros(3), space=sp)
    b = _OtherDummyField(values=np.zeros(3), space=sp)
    with pytest.raises(TypeError, match="same-class partner"):
        _ = a + b


@pytest.mark.foundation
def test_cross_class_subtraction_rejected():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.zeros(3), space=sp)
    b = _OtherDummyField(values=np.zeros(3), space=sp)
    with pytest.raises(TypeError, match="same-class partner"):
        _ = a - b


@pytest.mark.foundation
def test_cross_class_addition_rejected_on_shared_space():
    """The class-identity gate is units-independent.

    Two subclasses sharing the *same* geometric space (hence, under
    View-G, the same would-be units) MUST still refuse to add. This
    pins the "same units gives PERMISSION but not MEANING" semantics —
    same units is necessary-but-not-sufficient; class identity is what
    carries the physical kind. (Units are a role-leaf class constant,
    not a space property — issues #205 / #207.)
    """
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.zeros(3), space=sp)
    b = _OtherDummyField(values=np.zeros(3), space=sp)
    with pytest.raises(TypeError, match="same-class partner"):
        _ = a + b


# ─────────────────────────────────────────────────────────────────────
# Cross-space rejection (within same class)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_cross_space_addition_rejected_on_shape_mismatch():
    sp_a = FunctionSpace(name="dummy", shape=(3,))
    sp_b = FunctionSpace(name="dummy", shape=(4,))
    a = _DummyField(values=np.zeros(3), space=sp_a)
    b = _DummyField(values=np.zeros(4), space=sp_b)
    with pytest.raises(ValueError, match="equal space"):
        _ = a + b


@pytest.mark.foundation
def test_cross_space_addition_rejected_on_name_mismatch():
    sp_a = FunctionSpace(name="dummy_a", shape=(3,))
    sp_b = FunctionSpace(name="dummy_b", shape=(3,))
    a = _DummyField(values=np.zeros(3), space=sp_a)
    b = _DummyField(values=np.zeros(3), space=sp_b)
    with pytest.raises(ValueError, match="equal space"):
        _ = a + b


# ─────────────────────────────────────────────────────────────────────
# Diagnostics: linf, l2, inner_product, copy
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_linf_returns_max_absolute_value():
    sp = FunctionSpace(name="dummy", shape=(4,))
    a = _DummyField(values=np.array([1.0, -3.0, 2.0, -0.5]), space=sp)
    assert a.linf == pytest.approx(3.0)


@pytest.mark.foundation
def test_l2_reads_space_metric_euclidean():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.array([1.0, 2.0, 2.0]), space=sp)
    # Euclidean ‖a‖₂ = √(1 + 4 + 4) = 3.
    assert a.l2 == pytest.approx(3.0)


@pytest.mark.foundation
def test_l2_reads_space_metric_weighted():
    w = np.array([1.0, 2.0, 3.0])
    sp = FunctionSpace(name="dummy", shape=(3,), inner_product_weights=w)
    a = _DummyField(values=np.array([1.0, 1.0, 1.0]), space=sp)
    # ‖a‖_w = √(1·1 + 2·1 + 3·1) = √6.
    assert a.l2 == pytest.approx(np.sqrt(6.0))


@pytest.mark.foundation
def test_inner_product_euclidean_default():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.array([1.0, 2.0, 3.0]), space=sp)
    b = _DummyField(values=np.array([4.0, 5.0, 6.0]), space=sp)
    # 1·4 + 2·5 + 3·6 = 32.
    assert a.inner_product(b) == pytest.approx(32.0)


@pytest.mark.foundation
def test_inner_product_rejects_cross_class():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.zeros(3), space=sp)
    b = _OtherDummyField(values=np.zeros(3), space=sp)
    with pytest.raises(TypeError, match="same-class partner"):
        _ = a.inner_product(b)


@pytest.mark.foundation
def test_copy_returns_independent_ndarray():
    sp = FunctionSpace(name="dummy", shape=(3,))
    a = _DummyField(values=np.array([1.0, 2.0, 3.0]), space=sp)
    b = a.copy()
    assert type(b) is _DummyField
    assert b.space is sp
    np.testing.assert_array_equal(a.values, b.values)
    # Mutating one must not affect the other. The ABC is frozen, but
    # the underlying ndarray is owned per instance.
    b.values[0] = 99.0
    assert a.values[0] == 1.0


# ─────────────────────────────────────────────────────────────────────
# Subclass with extra fields — replace() propagation
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_subclass_with_extra_field_propagates_through_addition():
    """Subclasses that add fields (e.g. ``tag``, future ``mesh``)
    MUST have those fields preserved across dunder arithmetic
    without manual override. ``dataclasses.replace`` handles this
    automatically when ``kw_only=True``.
    """
    sp = FunctionSpace(name="rich", shape=(3,))
    a = _RichField(values=np.array([1.0, 2.0, 3.0]), space=sp, tag="alpha")
    b = _RichField(values=np.array([10.0, 20.0, 30.0]), space=sp, tag="alpha")
    c = a + b
    assert type(c) is _RichField
    assert c.tag == "alpha"
    assert c.space is sp
    np.testing.assert_array_equal(c.values, [11.0, 22.0, 33.0])


@pytest.mark.foundation
def test_subclass_with_extra_field_propagates_through_scalar_mul():
    sp = FunctionSpace(name="rich", shape=(2,))
    a = _RichField(values=np.array([1.0, 2.0]), space=sp, tag="beta")
    c = 3.0 * a
    assert type(c) is _RichField
    assert c.tag == "beta"
    np.testing.assert_array_equal(c.values, [3.0, 6.0])


@pytest.mark.foundation
def test_subclass_with_extra_field_propagates_through_negation_and_copy():
    sp = FunctionSpace(name="rich", shape=(2,))
    a = _RichField(values=np.array([1.0, 2.0]), space=sp, tag="gamma")
    neg = -a
    cp = a.copy()
    assert type(neg) is _RichField and neg.tag == "gamma"
    assert type(cp) is _RichField and cp.tag == "gamma"
    np.testing.assert_array_equal(neg.values, [-1.0, -2.0])
    np.testing.assert_array_equal(cp.values, [1.0, 2.0])


@pytest.mark.foundation
def test_where_largest_locates_the_peak():
    """where_largest(k) returns the k index tuples of largest |values|,
    largest first — the per-entry magnitude map (promoted from the retired
    displacement diagnostics surface at campaign 1 CS3, 2026-08-19)."""
    rng = np.random.default_rng(7)
    values = rng.normal(size=(4, 3, 5))
    field = _DummyField(
        values=values, space=FunctionSpace(name="map", shape=(4, 3, 5)),
    )
    flat = np.abs(values)
    peak = tuple(int(i) for i in np.unravel_index(int(flat.argmax()), flat.shape))
    got = field.where_largest(1)
    if got[0] != peak:
        raise AssertionError(f"where_largest peak {got[0]} != argmax {peak}")
    top3 = field.where_largest(3)
    if len(top3) != 3:
        raise AssertionError("where_largest(3) did not return 3 indices")
    mags = [flat[t] for t in top3]
    if mags != sorted(mags, reverse=True):
        raise AssertionError(f"where_largest(3) not largest-first: {mags}")
