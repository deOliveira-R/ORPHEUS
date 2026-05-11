r"""Tests for :func:`orpheus.sn.boundary_realize.realize_recursively`.

The walker is structural (Wave 11 / C11.2): it traverses a
:class:`~orpheus.numerics.operator.LinearOperator` expression tree
whose leaves are :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
instances and whose internal nodes are recognised Wave-0 composers
(``ScaledOperator``, ``OperatorSum``, ``OperatorProduct``,
``TensorProductOperator``, ``SumOfTensorProductsOperator``). Pure
tree-walking — no math content. Tests pin:

* **L0 / foundation** — walker correctness:

  - Leaf delegation: a bare
    :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` realises
    via :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`
    (identity-of-result with direct realisation).
  - Composer preservation: every recognised Wave-0 composer node
    survives the recursion with the same type, only the leaves replaced.
  - Nested expressions: nested ``ScaledOperator``/``OperatorSum``
    trees are walked depth-first.
  - Unknown node raises :class:`BoundaryError` naming the offending
    type.

* **L1 / apply equivalence** — the realised tree's ``apply(psi)``
  output equals the explicit pointwise weighted sum of leaf
  realisations applied independently (the structurally-independent
  reference for the Wave-0 ``OperatorSum`` /``ScaledOperator``
  algebra). The composition surface IS Wave-0; we verify against the
  pointwise sum (which uses only ``numpy`` addition + scalar
  multiplication, both above the trusted-library line per
  ``algebra-of-record``).

V&V tag: ``@pytest.mark.l0`` / ``@pytest.mark.foundation`` for the
walker-structure tests; ``@pytest.mark.l1`` for the apply-equivalence
tests. The bridge between L0 and L1 is the single-leaf identity
(L0 test ``test_leaf_realisation_matches_direct_realize``).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundaryOperator,
    BoundaryError,
    SpecularBoundaryOperator,
    VacuumBoundaryOperator,
    WhiteBoundaryOperator,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    OperatorProduct,
    OperatorSum,
    PermutationOperator,
    ScaledOperator,
)
from orpheus.sn.angular_operator import AngularAverageOperator
from orpheus.sn.boundary_realize import realize_recursively
from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace
from orpheus.sn.quadrature import GaussLegendre1D, LebedevSphere, LevelSymmetricSN


# ─────────────────────────────────────────────────────────────────────
# 1. Leaf delegation
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_leaf_realisation_matches_direct_realize() -> None:
    """A bare :class:`BoundaryTraceLaw` leaf realises via the SN realiser.

    The walker's leaf case MUST be a no-op wrapper around
    :meth:`SNBoundaryRealizer.realize` — anything else would break the
    L0 → L1 bridge that single-leaf equivalence carries.
    """
    quad = GaussLegendre1D.create(8)
    ms = SNMethodSpace.minimal(quad)
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)

    walker_op = realize_recursively(spec, ms)
    direct_op = SNBoundaryRealizer().realize(spec, ms)

    # The walker returns the realiser's output verbatim at a leaf —
    # for the α=1 fast path both are the bare PermutationOperator,
    # so test by apply-identity on a random psi (the realiser does
    # NOT memoize, so identity-of-object is not guaranteed).
    rng = np.random.default_rng(0)
    psi = rng.standard_normal((quad.N, 3))
    np.testing.assert_array_equal(walker_op.apply(psi), direct_op.apply(psi))


@pytest.mark.foundation
def test_leaf_vacuum_realisation_matches_direct_realize() -> None:
    """Vacuum (which requires ``inflow_indices``) realises via the leaf path."""
    quad = LebedevSphere.create(17)
    inflow_indices = np.flatnonzero(quad.mu_x > 0)
    ms = SNMethodSpace(
        quadrature=quad, face="xmin", inflow_indices=inflow_indices,
    )

    walker_op = realize_recursively(VacuumBoundaryOperator(), ms)
    direct_op = SNBoundaryRealizer().realize(VacuumBoundaryOperator(), ms)

    rng = np.random.default_rng(1)
    psi = rng.standard_normal((quad.N, 4))
    np.testing.assert_array_equal(walker_op.apply(psi), direct_op.apply(psi))


# ─────────────────────────────────────────────────────────────────────
# 2. ScaledOperator composer
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_scaled_operator_around_leaf_preserves_composer() -> None:
    r"""``0.5 * ReflectiveBoundary("x")`` walks to ``ScaledOperator(0.5, P)``.

    The walker's ``ScaledOperator`` branch MUST preserve the scalar
    and wrap the realised inner operator. At α_BC=1 the inner realisation
    is the bare :class:`PermutationOperator`.
    """
    quad = LevelSymmetricSN.create(4)
    ms = SNMethodSpace.minimal(quad)
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
    law = 0.5 * spec

    op = realize_recursively(law, ms)
    assert isinstance(op, ScaledOperator)
    assert op.scalar == 0.5
    # Inner is the bare PermutationOperator (α=1 fast path in the
    # realiser); preserved verbatim by the walker.
    assert isinstance(op.op, PermutationOperator)


@pytest.mark.foundation
def test_scaled_operator_apply_matches_leaf_realisation_scaled() -> None:
    r"""``(c * leaf).apply(psi) == c * leaf_realised.apply(psi)`` to nulp=4."""
    quad = GaussLegendre1D.create(8)
    ms = SNMethodSpace.minimal(quad)
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)

    realised_leaf = SNBoundaryRealizer().realize(spec, ms)
    composed = realize_recursively(0.3 * spec, ms)

    rng = np.random.default_rng(2)
    psi = rng.standard_normal((quad.N, 2))
    np.testing.assert_array_almost_equal_nulp(
        composed.apply(psi), 0.3 * realised_leaf.apply(psi), nulp=4,
    )


# ─────────────────────────────────────────────────────────────────────
# 3. OperatorSum composer (the Wave-11 replacement for MixedBoundaryOperator)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
def test_operator_sum_of_scaled_leaves_realises_to_operator_sum() -> None:
    r"""Standard Marshak mixed: ``0.3 * Reflective + 0.7 * White``.

    The walker MUST recurse into the binary ``OperatorSum`` operands
    independently and reassemble via the same composer; each operand
    is itself a ``ScaledOperator`` around a leaf, so the realised
    output is ``OperatorSum(ScaledOperator(0.3, perm),
    ScaledOperator(0.7, angular_avg))``.
    """
    quad = LebedevSphere.create(17)
    ms = SNMethodSpace.minimal(quad)
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
    white = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
    law = 0.3 * spec + 0.7 * white

    op = realize_recursively(law, ms)
    assert isinstance(op, OperatorSum)
    assert isinstance(op.a, ScaledOperator)
    assert isinstance(op.b, ScaledOperator)
    assert op.a.scalar == 0.3
    assert op.b.scalar == 0.7
    # The inner realisations are bare primitives (α=1 fast paths).
    assert isinstance(op.a.op, PermutationOperator)
    assert isinstance(op.b.op, AngularAverageOperator)


@pytest.mark.l1
def test_operator_sum_apply_matches_explicit_weighted_sum() -> None:
    r"""``walker(0.3*spec + 0.7*white).apply(psi)`` matches the
    pointwise weighted sum of leaf realisations applied independently.

    The pointwise sum is the structurally-independent reference (its
    composition uses only ``numpy`` addition + scalar multiplication,
    above the trusted-library line per ``algebra-of-record``). nulp=4
    margin covers the binary-reduction FP order vs the pointwise
    reduction.
    """
    quad = LebedevSphere.create(17)
    ms = SNMethodSpace.minimal(quad)
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
    white = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)

    spec_realised = SNBoundaryRealizer().realize(spec, ms)
    white_realised = SNBoundaryRealizer().realize(white, ms)

    composed = realize_recursively(0.3 * spec + 0.7 * white, ms)

    rng = np.random.default_rng(3)
    psi = rng.uniform(0.0, 2.0, size=(quad.N, 5, 3))
    expected = 0.3 * spec_realised.apply(psi) + 0.7 * white_realised.apply(psi)

    np.testing.assert_array_almost_equal_nulp(
        composed.apply(psi), expected, nulp=4,
    )


# ─────────────────────────────────────────────────────────────────────
# 4. OperatorProduct composer
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_operator_product_of_albedo_and_identity_realises() -> None:
    r"""``AlbedoBoundary(0.5) @ AlbedoBoundary(0.4)`` realises to a
    binary ``OperatorProduct``.

    The physical interpretation of "compose two pure-attenuator BCs"
    is a successive-attenuation scenario (rare in production transport
    but algebraically valid); the test pins the WALKER's structural
    behaviour, NOT a physics claim. The walker MUST recurse into both
    operands and reassemble via :class:`OperatorProduct`.
    """
    quad = GaussLegendre1D.create(4)
    ms = SNMethodSpace.minimal(quad)
    a = AlbedoBoundaryOperator(albedo=0.5)
    b = AlbedoBoundaryOperator(albedo=0.4)
    law = a @ b

    op = realize_recursively(law, ms)
    assert isinstance(op, OperatorProduct)
    # apply equivalence: (A @ B).apply(psi) = 0.5 * (0.4 * psi) = 0.2 * psi.
    rng = np.random.default_rng(4)
    psi = rng.standard_normal((quad.N, 2))
    np.testing.assert_allclose(op.apply(psi), 0.2 * psi, rtol=1e-15)


# ─────────────────────────────────────────────────────────────────────
# 5. Unknown node raises BoundaryError
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_unknown_node_type_raises_boundary_error() -> None:
    """A node neither :class:`BoundaryTraceLaw` nor recognised composer
    raises :class:`BoundaryError` naming the offending type."""
    quad = GaussLegendre1D.create(4)
    ms = SNMethodSpace.minimal(quad)

    class _NotAComposer:  # noqa: D401 — ad-hoc test stand-in
        pass

    with pytest.raises(BoundaryError) as excinfo:
        realize_recursively(_NotAComposer(), ms)
    assert excinfo.value.law == "_NotAComposer"


@pytest.mark.foundation
def test_numpy_array_as_leaf_raises_boundary_error() -> None:
    """A plain ndarray is neither a BoundaryTraceLaw nor a Wave-0
    composer — the walker MUST raise."""
    quad = GaussLegendre1D.create(4)
    ms = SNMethodSpace.minimal(quad)

    with pytest.raises(BoundaryError):
        realize_recursively(np.eye(3), ms)


# ─────────────────────────────────────────────────────────────────────
# 6. Nested composition
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_nested_scaled_around_sum_walks_depth_first() -> None:
    r"""``0.5 * (0.3 * spec + 0.7 * white)`` walks depth-first.

    The walker MUST recurse into the outer ``ScaledOperator``'s inner
    operand (an ``OperatorSum``), then into each summand
    (``ScaledOperator``-around-leaf), preserving every composer
    layer.
    """
    quad = LebedevSphere.create(17)
    ms = SNMethodSpace.minimal(quad)
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
    white = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
    law = 0.5 * (0.3 * spec + 0.7 * white)

    op = realize_recursively(law, ms)
    assert isinstance(op, ScaledOperator)
    assert op.scalar == 0.5
    inner_sum = op.op
    assert isinstance(inner_sum, OperatorSum)
    assert isinstance(inner_sum.a, ScaledOperator)
    assert isinstance(inner_sum.b, ScaledOperator)
    assert inner_sum.a.scalar == 0.3
    assert inner_sum.b.scalar == 0.7


@pytest.mark.l1
def test_nested_apply_matches_pointwise_distribution() -> None:
    r"""``walker(0.5 * (0.3 * spec + 0.7 * white)).apply(psi)`` matches
    ``0.5 * (0.3 * spec_realised.apply(psi) + 0.7 * white_realised.apply(psi))``.

    Verifies the walker preserves the scalar-distributes-over-sum
    algebraic identity through depth-first recursion. nulp=4 covers
    the binary-reduction FP order.
    """
    quad = LebedevSphere.create(17)
    ms = SNMethodSpace.minimal(quad)
    spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
    white = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)

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
