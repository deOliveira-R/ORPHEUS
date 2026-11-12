r"""#707 Step A1, the member seam's laws on the numerics spaces: a space mints its members, and its metric names its kernel.

Two verbs landed with A1 (``.claude/plans/true_ends.md``, "Compaction point 1", sub-step A1):

* ``FunctionSpace.member(c)``, the inverse of the coordinate chart ``coordinates_of``. The class follows the space
  (ruling R5): a primal leaf mints a bare :class:`~orpheus.numerics.space.Member`, a dual leaf a
  :class:`~orpheus.numerics.space.Covector`.
* ``FunctionSpace.kernel_basis()`` (through ``HilbertMetric.kernel_basis(shape)``), an orthonormal basis of
  :math:`\ker G` in the flat C-order chart, the data the #705 kernel-condition check applies an operator to.

The laws, one row family each (THEOREM rows unless marked; ``foundation``: software invariants of the space algebra,
no theory-page label):

* M: ``coordinates_of(member(c)) == c`` with the space's shape, ``member(c).space is V``, the class R5 names; and the
  other leg, ``member(coordinates_of(x))`` equal to ``x`` in coordinates and class for an ``x`` built by its own
  constructor (never by ``member``, so the two legs are independent).
* K: the columns lie in :math:`\ker G` (``G`` read column by column off the space's own ``apply_metric`` on unit
  vectors), they are orthonormal (:math:`K^\top K = I`), and their count is the kernel dimension counted BY HAND from
  the fixture's construction (never from ``V.rank``, which reads the same stored data ``kernel_basis`` reads); a
  full-rank space has zero columns.
* KC: a direct sum's basis is block-wise, in the flat chart's block order: each column is supported on exactly one
  block, and there it lies in that block's kernel.

The production builders' composites (SN, diffusion) are posed in
``tests/gates/cross_method/test_member_seam_builders.py``: a numerics gate may not import them (the test-layer gate).

Tolerances: ``64 ε ‖G‖₂ ‖k‖`` for ``G k = 0`` and ``64 ε`` for ``KᵀK = I`` (an ``eigh`` of a 3 × 3 projector and
``kron`` products of at most three factors; the unit-column diagonal cases are exact and pass at any tolerance).
``-O``-safe: every check is ``np.testing`` or ``pytest.fail``.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import pytest

from orpheus.numerics.axis import Axis
from orpheus.numerics.coupled_system import CoupledSpace
from orpheus.geometry.manifold import RealSpace
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.metric import DenseMetric, DiagonalMetric, FactoredMetric
from orpheus.numerics.space import Covector, FunctionSpace, Member
from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from tests.gates.numerics._axis_marginal_lift_support import zero_weight_nodal

pytestmark = pytest.mark.foundation

_EPS = float(np.finfo(float).eps)
_EXACT = 64.0 * _EPS

_NM1 = "tests/gates/numerics/test_true_ends_member_seam.py::test_nm1_a_space_mints_its_members"
_RANK = "tests/gates/numerics/test_quotient_metric.py::test_B1_metric_rank_is_the_hand_count"

#: The rank-2 Gram on ℝ³ with exact integer entries, ``G = B Bᵀ``; its kernel is spanned by (1, −1, −1).
_B = np.array([[1.0, 1.0], [0.0, 1.0], [1.0, 0.0]])


def _axis(label: str, weights) -> Axis:
    w = np.asarray(weights, dtype=float)
    return Axis.of_measure(DiscreteMeasure(nodes=np.arange(w.size, dtype=float), weights=w, support=RealSpace(1)), label)


def _two_axis() -> FunctionSpace:
    """NM-1's fixture: weights 0.5, 2, 4 on one axis, 1 and 3 on the other (a stray metric factor is visible)."""
    return FunctionSpace.of_axes(_axis("a", [0.5, 2.0, 4.0]), _axis("b", [1.0, 3.0]))


def _diag(name: str, weights) -> FunctionSpace:
    return FunctionSpace(name, (len(weights),), inner_product_weights=np.asarray(weights, dtype=float))


def _dense_singular(name: str) -> FunctionSpace:
    return FunctionSpace(name, (3,), metric=DenseMetric(_B @ _B.T))


def _dense_near_cutoff(name: str) -> FunctionSpace:
    r"""A 6×6 Gram with spectrum (1, 0.5, 0.2, 0.1, 1e-8, 0) under a fixed rotation: one kernel direction beside a
    live eigenvalue far below the others (the A1 qa review's probe, 2026-11-12). An eigen-solve of the projector
    :math:`I − G^{+}G` loses digits as :math:`\varepsilon/\lambda_{\min}` here (‖GK‖/‖G‖ ≈ 1e-10); the eigenvectors
    of :math:`G` stay at rounding."""
    q, _ = np.linalg.qr(np.random.default_rng(29).standard_normal((6, 6)))
    gram = q @ np.diag([1.0, 0.5, 0.2, 0.1, 1e-8, 0.0]) @ q.T
    return FunctionSpace(name, (6,), metric=DenseMetric((gram + gram.T) / 2.0))


def _dense_full(name: str) -> FunctionSpace:
    return FunctionSpace(name, (3,), metric=DenseMetric(np.array([[2.0, 0.5, 0.0], [0.5, 1.0, 0.25], [0.0, 0.25, 4.0]])))


def _mixed_factored() -> FunctionSpace:
    r"""The afterthought: a degenerate diagonal block (weights 1, 0, 2) tensor the singular dense block ``B Bᵀ``.

    Both terms of :math:`\ker(G_1 \otimes G_2) = \ker G_1 \otimes V_2 + \operatorname{range} G_1 \otimes \ker G_2`
    are active (3 + 2 columns of 9), which no other fixture in the tree activates (``zero_weight_nodal`` is
    diagonal ⊗ diagonal, where the second term is a set of unit vectors like the first).
    """
    return FunctionSpace(
        "kb_mixed", (3, 3),
        metric=FactoredMetric((((3,), DiagonalMetric(np.array([1.0, 0.0, 2.0]))), ((3,), DenseMetric(_B @ _B.T)))),
    )


# ═══════════════════════════════════════════════════════════════════════
# M — member is the inverse of the chart, and the class follows the space
# ═══════════════════════════════════════════════════════════════════════

#: name -> (space builder, the class R5 names for its members).
_MEMBER_SPACES: dict[str, tuple[Callable[[], FunctionSpace], type]] = {
    "leaf diagonal (two axes)": (_two_axis, Member),
    "leaf diagonal, zero weight": (lambda: _diag("m_z", [1.0, 0.0, 2.0]), Member),
    "leaf dense": (lambda: _dense_full("m_d"), Member),
    "leaf factored zero-weight": (zero_weight_nodal, Member),
    "dual of leaf diagonal (two axes)": (lambda: _two_axis().dual(), Covector),
    "dual of leaf dense": (lambda: _dense_full("m_dd").dual(), Covector),
    "tensor product": (lambda: _diag("m_t1", [1.0, 0.0, 2.0]) * _diag("m_t2", [4.0, 0.5]), Member),
}


@pytest.mark.rests_on(_NM1)
@pytest.mark.parametrize("case", list(_MEMBER_SPACES))
def test_m1_coordinates_of_member_is_the_identity(case: str) -> None:
    """``coordinates_of(member(c)) == c`` (bit-exact: a reshape and a float cast), in the space's own shape; the member
    lives on the space (``is``) and is of the class R5 names. Fed both a flat and a shaped ``c``.
    First red: ``member`` minting on another space, or ``DualSpace.member`` returning a bare ``Member``."""
    builder, cls = _MEMBER_SPACES[case]
    V = builder()
    c = np.random.default_rng(7).standard_normal(V.n_coordinates)
    for given in (c, c.reshape(V.shape)):
        x = V.member(given)
        if type(x) is not cls:
            pytest.fail(f"[{case}] member class {type(x).__name__}, R5 names {cls.__name__}")
        if x.space is not V:
            pytest.fail(f"[{case}] member(c).space is not V")
        coords = np.asarray(V.coordinates_of(x))
        if coords.shape != V.shape:
            pytest.fail(f"[{case}] coordinates_of shape {coords.shape}, the space's {V.shape}")
        np.testing.assert_array_equal(coords.ravel(), c, err_msg=f"[{case}] coordinates_of(member(c)) != c")


@pytest.mark.rests_on(_NM1)
@pytest.mark.parametrize("case", list(_MEMBER_SPACES))
def test_m2_member_of_coordinates_is_the_identity(case: str) -> None:
    """The other leg: ``x`` built by its class's own constructor (``Member(values=…, space=V)``), then
    ``member(coordinates_of(x))`` equals ``x`` in coordinates (bit-exact), in class, and lives on V.
    First red: ``member`` permuting or scaling the coordinates it is handed."""
    builder, cls = _MEMBER_SPACES[case]
    V = builder()
    values = np.random.default_rng(11).standard_normal(V.shape)
    x = cls(values=values, space=V)
    y = V.member(V.coordinates_of(x))
    if type(y) is not type(x):
        pytest.fail(f"[{case}] member(coordinates_of(x)) is a {type(y).__name__}, x a {type(x).__name__}")
    if y.space is not V:
        pytest.fail(f"[{case}] member(coordinates_of(x)).space is not V")
    np.testing.assert_array_equal(np.asarray(y.values), values, err_msg=f"[{case}] member(coordinates_of(x)) != x")


# ═══════════════════════════════════════════════════════════════════════
# K — kernel_basis spans ker G, orthonormally, with the hand count
# ═══════════════════════════════════════════════════════════════════════

def _two_axis_block_factored() -> FunctionSpace:
    """A factored metric whose first block is a leading-aligned diagonal over a two-axis ``(3, 2)`` block (the
    elegance review's probe 4, 2026-11-12), tensored with a full-rank diagonal on ``(2,)``."""
    block = DiagonalMetric(np.array([1.0, 0.0, 2.0]))
    metric = FactoredMetric((((3, 2), block), ((2,), DiagonalMetric(np.array([1.0, 3.0])))))
    return FunctionSpace("kb_fb", (3, 2, 2), metric=metric)


#: name -> (space builder, the kernel dimension counted by hand from the construction).
_KB_SPACES: dict[str, tuple[Callable[[], FunctionSpace], int]] = {
    "diag one zero": (lambda: _diag("kb_d1", [1.0, 0.0, 2.0]), 1),
    "diag two zeros": (lambda: _diag("kb_d2", [0.0, 1.0, 0.0]), 2),
    "dense rank 2": (lambda: _dense_singular("kb_x"), 1),
    "dual of diag one zero": (lambda: _diag("kb_dd", [1.0, 0.0, 2.0]).dual(), 1),
    "dual of dense rank 2": (lambda: _dense_singular("kb_xd").dual(), 1),
    # First red: the kernel read off the projector I − G⁺G, whose columns leave ker G by ≈ 1e-10 ≫ 64ε.
    "dense near the cutoff": (lambda: _dense_near_cutoff("kb_nc"), 1),
    # 3 x 4, one zero angular weight: 1 * 4 kernel columns.
    "factored zero-weight nodal": (zero_weight_nodal, 4),
    # 9 = 3 x 3: rank 2 * 2 = 4, so 5 kernel columns (3 from ker G1 ⊗ V2, 2 from range G1 ⊗ ker G2).
    "factored diagonal ⊗ singular dense": (_mixed_factored, 5),
    # 6 = 3 x 2, the first factor one zero weight, the second full rank: 1 * 2.
    "tensor product diag ⊗ diag": (lambda: _diag("kb_t1", [1.0, 0.0, 2.0]) * _diag("kb_t2", [4.0, 0.5]), 2),
    # 12 = (3 x 2) x 2, a diagonal block over TWO axes whose weights lead (axis 0): one zero weight leaves 2 x 2 = 4
    # kernel columns. First red: the block's weights broadcast numpy's way (trailing), which refuses (3,) on (3, 2).
    "factored two-axis diagonal block": (_two_axis_block_factored, 4),
    "full rank diag": (lambda: _diag("kb_f", [1.0, 2.0, 0.5]), 0),
    "full rank dense": (lambda: _dense_full("kb_fd"), 0),
    "full rank two axes": (_two_axis, 0),
    "Euclidean": (lambda: FunctionSpace("kb_e", (4,)), 0),
}


def _gram(V: FunctionSpace) -> np.ndarray:
    """``G`` column by column off the space's own ``apply_metric`` on the unit vectors of the flat chart."""
    n = V.n_coordinates
    return np.column_stack([np.asarray(V.apply_metric(e.reshape(V.shape)), dtype=float).ravel() for e in np.eye(n)])


def _check_kernel_columns(case: str, G: np.ndarray, K: np.ndarray) -> None:
    norm_G = np.linalg.norm(G, 2)
    for j in range(K.shape[1]):
        if np.linalg.norm(G @ K[:, j]) > _EXACT * norm_G * np.linalg.norm(K[:, j]):
            pytest.fail(f"[{case}] column {j} is not in ker G: ‖G k‖ = {np.linalg.norm(G @ K[:, j]):.2e}")
    gram = K.T @ K
    if np.abs(gram - np.eye(K.shape[1])).max(initial=0.0) > _EXACT:
        pytest.fail(f"[{case}] the kernel basis is not orthonormal: max|KᵀK − I| = {np.abs(gram - np.eye(K.shape[1])).max():.2e}")


@pytest.mark.rests_on(_RANK)
@pytest.mark.parametrize("case", list(_KB_SPACES))
def test_k1_the_kernel_basis_spans_the_kernel(case: str) -> None:
    """The four laws on a leaf, a dual, a factored and a tensor-product space: shape ``(n, hand count)``; every column
    in ``ker G``; ``KᵀK = I``; a full-rank space has zero columns. In-kernel, orthonormal and of the kernel's dimension
    is a basis of the kernel. The count is also ``n − rank`` (the two production counts must agree).
    First reds: a complemented zero mask (columns in the range), the dense block's live and kernel bases swapped, the
    factored product's second term dropped (3 columns for 5)."""
    builder, hand = _KB_SPACES[case]
    V = builder()
    n = V.n_coordinates
    K = np.asarray(V.kernel_basis(), dtype=float)
    if K.shape != (n, hand):
        pytest.fail(f"[{case}] kernel basis shape {K.shape}, the hand count gives {(n, hand)}")
    if hand != n - V.rank:
        pytest.fail(f"[{case}] the hand count {hand} disagrees with n − rank = {n - V.rank}")
    _check_kernel_columns(case, _gram(V), K)


# ═══════════════════════════════════════════════════════════════════════
# KC — a direct sum's kernel basis is block-wise, in the flat chart's order
# ═══════════════════════════════════════════════════════════════════════


def _composites() -> dict[str, tuple[FunctionSpace, tuple[FunctionSpace, ...], tuple[int, ...]]]:
    """name -> (composite, its blocks in flat-chart order, the hand kernel count per block). Block sizes and kernels
    differ (3 with kernel e₁, 3 dense with kernel (1, −1, −1), 12 with 4), so a reversed or shifted order is visible."""
    diag, dense, nodal = _diag("kc_v", [1.0, 0.0, 2.0]), _dense_singular("kc_d"), zero_weight_nodal()
    return {
        "FullFieldSpace (diag ⊕ dense)": (FullFieldSpace.from_blocks(diag, dense), (diag, dense), (1, 1)),
        "FullFieldSpace (nodal ⊕ diag)": (FullFieldSpace.from_blocks(nodal, diag), (nodal, diag), (4, 1)),
        "CoupledSpace (diag ⊕ dense ⊕ nodal)": (CoupledSpace.from_systems((diag, dense, nodal)), (diag, dense, nodal), (1, 1, 4)),
        "CoupledSpace (full rank ⊕ diag)": (
            CoupledSpace.from_systems((_dense_full("kc_f"), diag)), (_dense_full("kc_f2"), diag), (0, 1),
        ),
    }


@pytest.mark.rests_on(_RANK)
@pytest.mark.parametrize("case", list(_composites()))
def test_kc_a_direct_sum_kernel_basis_is_block_wise(case: str) -> None:
    """Shape ``(Σ n_i, Σ hand_i)``; ``KᵀK = I``; each column is supported on exactly one block (its other slices are
    exactly zero), lies in that block's kernel (the block's own ``G``), and block ``i`` holds ``hand_i`` columns.
    The composite's ``G`` is read per block because a composite's ``apply_metric`` takes a composite field.
    First red: the blocks' bases stacked in the reverse order."""
    V, blocks, hand = _composites()[case]
    K = np.asarray(V.kernel_basis(), dtype=float)
    sizes = [b.n_coordinates for b in blocks]
    if K.shape != (sum(sizes), sum(hand)):
        pytest.fail(f"[{case}] kernel basis shape {K.shape}, the hand count gives {(sum(sizes), sum(hand))}")
    edges = np.cumsum([0, *sizes])
    grams = [_gram(b) for b in blocks]
    per_block = [0] * len(blocks)
    for j in range(K.shape[1]):
        support = [i for i in range(len(blocks)) if np.any(K[edges[i]:edges[i + 1], j] != 0.0)]
        if len(support) != 1:
            pytest.fail(f"[{case}] column {j} is supported on blocks {support}, not on exactly one")
        i = support[0]
        per_block[i] += 1
        piece = K[edges[i]:edges[i + 1], j]
        if np.linalg.norm(grams[i] @ piece) > _EXACT * np.linalg.norm(grams[i], 2) * np.linalg.norm(piece):
            pytest.fail(f"[{case}] column {j} is not in block {i}'s kernel")
    if tuple(per_block) != hand:
        pytest.fail(f"[{case}] columns per block {tuple(per_block)}, the hand count gives {hand}")
    if np.abs(K.T @ K - np.eye(K.shape[1])).max(initial=0.0) > _EXACT:
        pytest.fail(f"[{case}] the kernel basis is not orthonormal")
