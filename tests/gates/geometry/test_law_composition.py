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
  :func:`~orpheus.transport.boundary_realization.realize_recursively` (the
  method-blind walker; its home since layer-DAG WP-E1, in geometry/boundary/
  from #290 P7b until then, with the leaf
  realizer — here :class:`SNBoundaryRealizer` — passed explicitly)
  transforms a descriptor tree into an operator tree:
  - The realised operator is the law tree's image in the operator algebra:
    LawSum and LawScaled realise through ``+`` and ``α·`` to one formal linear
    combination over the realised leaves (an ``OperatorSum``, flattened and
    with the scalars distributed, since S2), a leaf to a 1-arg LinearOperator.
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

import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    LawScaled,
    LawSum,
    ReflectiveBoundary,
    SpecularReturn,
    VacuumInflow,
    WhiteBoundary,
)

#: The specular leaf of every tree below is the perfect specular wall
#: ``AlbedoBoundary(1.0, SpecularReturn("x"))``, which realizes to the mirror's
#: bare tensor product. It was the mirror ``ReflectiveBoundary("x")`` itself
#: until the reflective cleanup (2026-10-01), which made a deck law (a
#: symmetry of the domain) refuse to be scaled or mixed; that refusal is
#: ``tests/gates/geometry/test_deck_laws_do_not_compose.py``.
#:
#: ⚠ The mixed-law gates below compose the specular wall with **white**, and
#: both are narrowed (B3.2 and B3.4a respectively), so the sum is well-typed and the
#: distributivity claim is stateable again — the two ``xfail(strict=True)``
#: markers that stood here were deleted when B3.4a landed, exactly as their own
#: reason text instructed.
#:
#: ⛔ **CORRECTED 2026-08-14 — the blocker below is DISCHARGED.** This read
#: "they must NOT be re-posed over **albedo** or **periodic**; those two are
#: still full-face endomorphisms (B3.4b / B3.4c) … a gate here over
#: ``0.3 * spec + 0.7 * albedo`` would be green and worthless until B3.4b
#: lands." **B3.4b (``943b37c1``) and B3.4c (``79b5affe``) both landed
#: 2026-08-01**, and `[M]` both laws narrow at HEAD: albedo through
#: ``_outflow_restriction(method_space, "albedo")``
#: (``sn/boundary/realizer.py:970``), periodic to :math:`\Gamma_+(\text{partner})
#: \to \Gamma_-(\text{face})` (``:1049``).
#:
#: The hazard the old note described was real and is worth keeping in view:
#: because ``|Γ₊| == |Γ₋|`` on every quadrature × face in the tree, a sum
#: mixing an un-narrowed leaf with a narrowed one does not RAISE — it runs and
#: is silently wrong (vv Mode 12). That is why the gate had to wait. It no
#: longer does: **a mixed ``0.3 * spec + 0.7 * albedo`` gate is now stateable
#: and is not yet written.**


# ═════════════════════════════════════════════════════════════════════
# 1. Descriptor-tree algebra closure
# ═════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_law_plus_law_returns_lawsum() -> None:
    """``law + law`` returns a :class:`LawSum`."""
    spec = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    tree = spec + white
    assert isinstance(tree, LawSum)
    assert tree.a is spec
    assert tree.b is white


@pytest.mark.foundation
def test_scalar_times_law_returns_lawscaled() -> None:
    """``α * law`` returns a :class:`LawScaled` (and ``law * α`` too)."""
    spec = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
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
    spec = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
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
    spec = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
    folded = 2.0 * (3.0 * spec)
    assert isinstance(folded, LawScaled)
    assert folded.scalar == 6.0
    # Inner is the leaf (not another LawScaled).
    assert isinstance(folded.inner, AlbedoBoundary)
    # Folds in either order:
    folded_other = (3.0 * spec) * 2.0
    assert isinstance(folded_other, LawScaled)
    assert folded_other.scalar == 6.0


@pytest.mark.foundation
def test_lawscaled_truediv_inverts_scalar() -> None:
    """``law / α = (1/α) * law``."""
    spec = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
    tree = spec / 4.0
    assert isinstance(tree, LawScaled)
    assert tree.scalar == 0.25
    assert tree.inner is spec


@pytest.mark.foundation
def test_lawsum_minus_law_uses_minus_one_scaled() -> None:
    """``a - b`` rewrites as ``LawSum(a, LawScaled(-1, b))``."""
    spec = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
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
    spec = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
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
    a = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
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
    spec = AlbedoBoundary(1.0, SpecularReturn(axis="x"))
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














# ═════════════════════════════════════════════════════════════════════
# 4. Type-error guard on unknown node types
# ═════════════════════════════════════════════════════════════════════



# ─────────────────────────────────────────────────────────────────────
# RETIRED at wave 2 S2b (retirement-audit D.14 demotion check): the eight ``realize_recursively`` rows (leaf dispatch, scaled wrap,
# sum, pointwise weighted sum, depth-first walk, nested distributive form, the unknown-node and ndarray-leaf TypeErrors) pinned the
# tree walker ``realize_recursively`` against ``SNBoundaryRealizer``. The walker and the realizer class are deleted (the carve
# retires the descriptor tree: a weighted mixture of walls is ONE ``WallCondition`` whose network is a sum of chains, realized by
# the table to a ``RealizedSum``). Their claims: the leaf dispatch and the weighted sum are
# ``tests/gates/sn/operators/test_sn_boundary_realizer.py::TestRealizeTheWeightedWalls`` (the sum of the images in the canonical
# order, weight 1 the bare kernel and a partial weight the scaled one); distributivity and nesting are the algebra of
# ``WallCondition`` itself, ``tests/gates/transport/test_wall_condition_algebra.py``; the TypeError rows have no subject (a
# declaration is a typed value, not a duck-typed tree). Section 3 of this file is gone with them; the descriptor tree rows above
# (``LawSum``/``LawScaled`` construction) stay until their classes retire at WF2.4.

