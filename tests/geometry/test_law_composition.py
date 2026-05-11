r"""Tests for the descriptor-tree composition algebra (Issue #186 / B3 + β2).

:class:`~orpheus.geometry.boundary._composition.LawSum` and
:class:`~orpheus.geometry.boundary._composition.LawScaled` form a
closed algebra over
``BoundaryTraceLaw | LawSum | LawScaled``. The tree is a **pure
descriptor** — no ``apply`` method; calling it before realization is
a static-type error.

These tests pin:

* **L0 / foundation** — algebra closure:
  - ``law + law`` returns :class:`LawSum`.
  - ``α * law`` returns :class:`LawScaled`.
  - ``LawScaled * β`` folds constants: ``β * (α * x) = (αβ) * x``.
  - ``LawSum / α``, ``-LawSum`` produce the expected nodes.
  - Descriptor trees carry NO ``apply`` method (the static-type
    invariant).

* **L1 / realisation** —
  :func:`~orpheus.sn.boundary_realize.realize_recursively`
  transforms a descriptor tree into an operator tree:
  - The output type family matches the input (LawSum → OperatorSum,
    LawScaled → ScaledOperator, leaf → realized 1-arg LinearOperator).
  - The realised tree's ``apply(psi)`` matches the explicit pointwise
    weighted sum of leaf realisations applied independently.
  - Nested composition is walked depth-first.
  - Unknown node types raise :class:`TypeError`.

The structurally-independent reference for the apply-equivalence
tests is the explicit pointwise weighted sum, which uses only
``numpy`` addition and scalar multiplication (above the
trusted-library line per ``algebra-of-record``).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    LawScaled,
    LawSum,
    ReflectiveBoundary,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.numerics.operator import OperatorSum, ScaledOperator
from orpheus.sn.boundary_realize import realize_recursively
from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace
from orpheus.sn.quadrature import GaussLegendre1D, LebedevSphere


# ═════════════════════════════════════════════════════════════════════
# 1. Descriptor-tree algebra closure
# ═════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_law_plus_law_returns_lawsum() -> None:
    """``law + law`` returns a :class:`LawSum`."""
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    tree = spec + white
    assert isinstance(tree, LawSum)
    assert tree.a is spec
    assert tree.b is white


@pytest.mark.foundation
def test_scalar_times_law_returns_lawscaled() -> None:
    """``α * law`` returns a :class:`LawScaled` (and ``law * α`` too)."""
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    left = 0.3 * spec
    right = spec * 0.3
    assert isinstance(left, LawScaled)
    assert isinstance(right, LawScaled)
    assert left.scalar == 0.3
    assert left.inner is spec
    assert right.scalar == 0.3
    assert right.inner is spec


@pytest.mark.foundation
def test_marshak_form_builds_law_sum_of_law_scaled() -> None:
    """``0.3 * spec + 0.7 * white`` is a :class:`LawSum` whose operands
    are :class:`LawScaled` wrappers around leaves."""
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    tree = 0.3 * spec + 0.7 * white
    assert isinstance(tree, LawSum)
    assert isinstance(tree.a, LawScaled)
    assert tree.a.scalar == 0.3
    assert tree.a.inner is spec
    assert isinstance(tree.b, LawScaled)
    assert tree.b.scalar == 0.7
    assert tree.b.inner is white


@pytest.mark.foundation
def test_lawscaled_constant_folding_collapses_chain() -> None:
    """``α * (β * law) = (α*β) * law`` — the chain never re-nests."""
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    folded = 2.0 * (3.0 * spec)
    assert isinstance(folded, LawScaled)
    assert folded.scalar == 6.0
    # Inner is the leaf (not another LawScaled).
    assert isinstance(folded.inner, ReflectiveBoundary)
    # Folds in either order:
    folded_other = (3.0 * spec) * 2.0
    assert isinstance(folded_other, LawScaled)
    assert folded_other.scalar == 6.0


@pytest.mark.foundation
def test_lawscaled_truediv_inverts_scalar() -> None:
    """``law / α = (1/α) * law``."""
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    tree = spec / 4.0
    assert isinstance(tree, LawScaled)
    assert tree.scalar == 0.25
    assert tree.inner is spec


@pytest.mark.foundation
def test_lawsum_minus_law_uses_minus_one_scaled() -> None:
    """``a - b`` rewrites as ``LawSum(a, LawScaled(-1, b))``."""
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    tree = spec - white
    assert isinstance(tree, LawSum)
    assert tree.a is spec
    assert isinstance(tree.b, LawScaled)
    assert tree.b.scalar == -1.0
    assert tree.b.inner is white


@pytest.mark.foundation
def test_neg_lawsum_wraps_as_minus_one_scaled() -> None:
    """``-sum_node`` wraps as ``LawScaled(-1, sum_node)``."""
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    s = spec + white
    neg = -s
    assert isinstance(neg, LawScaled)
    assert neg.scalar == -1.0
    assert neg.inner is s


@pytest.mark.foundation
def test_lawsum_plus_law_returns_lawsum() -> None:
    """``LawSum + law`` is closed: returns a deeper :class:`LawSum`.

    The tree is **not** flattened: ``(a + b) + c`` is
    ``LawSum(LawSum(a, b), c)``, distinct from ``LawSum(a, LawSum(b, c))``.
    """
    a = ReflectiveBoundary(axis="x", albedo=1.0)
    b = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    c = AlbedoBoundary(albedo=0.5)
    tree = (a + b) + c
    assert isinstance(tree, LawSum)
    # Left operand is the inner LawSum.
    assert isinstance(tree.a, LawSum)
    assert tree.a.a is a
    assert tree.a.b is b
    # Right operand is c.
    assert tree.b is c


# ═════════════════════════════════════════════════════════════════════
# 2. Descriptors carry NO apply method (the static-type invariant)
# ═════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_leaf_descriptor_has_no_apply() -> None:
    """A bare BC descriptor has no ``apply``. Issue #186 invariant."""
    assert not hasattr(ReflectiveBoundary(axis="x"), "apply")
    assert not hasattr(WhiteBoundary(axis="x", outward_sign=+1), "apply")
    assert not hasattr(VacuumInflow(), "apply")
    assert not hasattr(AlbedoBoundary(albedo=0.5), "apply")


@pytest.mark.foundation
def test_law_tree_has_no_apply() -> None:
    """LawSum / LawScaled compositions also have no ``apply``."""
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    sum_tree = spec + white
    scaled_tree = 0.5 * spec
    marshak = 0.3 * spec + 0.7 * white
    nested = 0.5 * (0.3 * spec + 0.7 * white)
    for tree in (sum_tree, scaled_tree, marshak, nested):
        assert not hasattr(tree, "apply"), (
            f"{type(tree).__name__} should not carry apply"
        )


# ═════════════════════════════════════════════════════════════════════
# 3. realize_recursively: descriptor tree → operator tree
# ═════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_realize_recursively_leaf_dispatches_to_sn_realizer() -> None:
    """A bare :class:`BoundaryTraceLaw` leaf realises via the SN realiser."""
    quad = GaussLegendre1D.create(8)
    ms = SNMethodSpace.minimal(quad)
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    walker_op = realize_recursively(spec, ms)
    direct_op = SNBoundaryRealizer().realize(spec, ms)
    rng = np.random.default_rng(0)
    psi = rng.standard_normal((quad.N, 3))
    np.testing.assert_array_equal(walker_op.apply(psi), direct_op.apply(psi))


@pytest.mark.foundation
def test_realize_recursively_lawscaled_wraps_in_scaled_operator() -> None:
    """``LawScaled(α, leaf)`` realises to ``ScaledOperator(α, realised_leaf)``."""
    quad = GaussLegendre1D.create(8)
    ms = SNMethodSpace.minimal(quad)
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    tree = 0.5 * spec
    op = realize_recursively(tree, ms)
    assert isinstance(op, ScaledOperator)
    assert op.scalar == 0.5
    # The inner is the realised leaf (PermutationOperator at α=1).
    rng = np.random.default_rng(1)
    psi = rng.standard_normal((quad.N, 2))
    direct = SNBoundaryRealizer().realize(spec, ms)
    np.testing.assert_array_almost_equal_nulp(
        op.apply(psi), 0.5 * direct.apply(psi), nulp=4,
    )


@pytest.mark.foundation
def test_realize_recursively_lawsum_returns_operator_sum() -> None:
    """``LawSum(a, b)`` realises to ``OperatorSum(realise(a), realise(b))``."""
    quad = LebedevSphere.create(17)
    ms = SNMethodSpace.minimal(quad)
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    tree = 0.3 * spec + 0.7 * white
    op = realize_recursively(tree, ms)
    assert isinstance(op, OperatorSum)
    assert isinstance(op.a, ScaledOperator)
    assert isinstance(op.b, ScaledOperator)
    assert op.a.scalar == 0.3
    assert op.b.scalar == 0.7


@pytest.mark.l1
def test_realize_recursively_apply_matches_pointwise_weighted_sum() -> None:
    """``walker(0.3*spec + 0.7*white).apply(psi)`` matches the pointwise
    weighted sum of leaf realisations applied independently.

    The pointwise sum is the structurally-independent reference (uses
    only numpy + and *, above the trusted-library line per
    ``algebra-of-record``). nulp=4 covers the binary-reduction FP
    order vs the pointwise reduction.
    """
    quad = LebedevSphere.create(17)
    ms = SNMethodSpace.minimal(quad)
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    spec_realised = SNBoundaryRealizer().realize(spec, ms)
    white_realised = SNBoundaryRealizer().realize(white, ms)
    composed = realize_recursively(0.3 * spec + 0.7 * white, ms)
    rng = np.random.default_rng(3)
    psi = rng.uniform(0.0, 2.0, size=(quad.N, 5, 3))
    expected = 0.3 * spec_realised.apply(psi) + 0.7 * white_realised.apply(psi)
    np.testing.assert_array_almost_equal_nulp(
        composed.apply(psi), expected, nulp=4,
    )


@pytest.mark.foundation
def test_realize_recursively_walks_nested_depth_first() -> None:
    """``0.5 * (0.3 * spec + 0.7 * white)`` walks depth-first.

    The outer ``LawScaled(0.5, ...)`` wraps a ``LawSum``; each summand
    is a ``LawScaled`` around a leaf. The realised tree mirrors the
    structure with Wave-0 composers.
    """
    quad = LebedevSphere.create(17)
    ms = SNMethodSpace.minimal(quad)
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    tree = 0.5 * (0.3 * spec + 0.7 * white)
    op = realize_recursively(tree, ms)
    assert isinstance(op, ScaledOperator)
    assert op.scalar == 0.5
    inner_sum = op.op
    assert isinstance(inner_sum, OperatorSum)
    assert isinstance(inner_sum.a, ScaledOperator)
    assert isinstance(inner_sum.b, ScaledOperator)
    assert inner_sum.a.scalar == 0.3
    assert inner_sum.b.scalar == 0.7


@pytest.mark.l1
def test_realize_recursively_nested_apply_matches_distributive_form() -> None:
    """``walker(0.5 * (0.3 * spec + 0.7 * white)).apply(psi)`` matches
    ``0.5 * (0.3 * spec_realised.apply(psi) + 0.7 * white_realised.apply(psi))``.

    Verifies the walker preserves the scalar-distributes-over-sum
    identity through depth-first recursion. ``nulp=4`` covers
    binary-reduction FP order.
    """
    quad = LebedevSphere.create(17)
    ms = SNMethodSpace.minimal(quad)
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    spec_realised = SNBoundaryRealizer().realize(spec, ms)
    white_realised = SNBoundaryRealizer().realize(white, ms)
    composed = realize_recursively(0.5 * (0.3 * spec + 0.7 * white), ms)
    rng = np.random.default_rng(5)
    psi = rng.uniform(0.0, 2.0, size=(quad.N, 4))
    expected = 0.5 * (
        0.3 * spec_realised.apply(psi)
        + 0.7 * white_realised.apply(psi)
    )
    np.testing.assert_array_almost_equal_nulp(
        composed.apply(psi), expected, nulp=4,
    )


# ═════════════════════════════════════════════════════════════════════
# 4. Type-error guard on unknown node types
# ═════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_realize_recursively_raises_type_error_on_unknown_node() -> None:
    """A node neither :class:`BoundaryTraceLaw` nor :class:`LawSum` /
    :class:`LawScaled` raises :class:`TypeError`."""
    quad = GaussLegendre1D.create(4)
    ms = SNMethodSpace.minimal(quad)

    class _NotALaw:
        pass

    with pytest.raises(TypeError, match="LawSum | LawScaled"):
        realize_recursively(_NotALaw(), ms)


@pytest.mark.foundation
def test_realize_recursively_raises_on_ndarray_leaf() -> None:
    """A plain ndarray is not a descriptor — walker raises
    :class:`TypeError`."""
    quad = GaussLegendre1D.create(4)
    ms = SNMethodSpace.minimal(quad)
    with pytest.raises(TypeError, match="LawSum | LawScaled"):
        realize_recursively(np.eye(3), ms)
