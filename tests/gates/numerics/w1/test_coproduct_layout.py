r"""G05p: the offsets of a layout are derived from its sizes and tile ``range(n)`` exactly (N05, spec 7.6 row G05p, ON05-6, law L13).

``CoproductLayout`` (``numerics/coproduct_layout.py``, the renamed ``FaceLayout``) takes ordered ``(key, shape)`` pairs; ``offset`` is ``init=False`` on ``LayoutBlock``, so a layout whose offsets disagree with
its sizes is unspellable and the "gap or overlap" refusal that policed it retires; this gate is the witness of the mechanism that replaced the refusal. The assertions are hand rules written from the
definitions: the expected offsets are the cumulative sums of the products of the stated shapes, and the tiling is read by counting, for each flat index, the blocks that cover it.

First red: import-only at the base tip (the module ``orpheus.numerics.coproduct_layout`` does not exist before the rename). Behavioural red at the tree of commit (2) (the renamed layout with explicit
offsets): ``LayoutBlock(offset=...)`` is accepted (and, after ruling V3, ``flat_size=...``, a derived property), and ``layout.offsets`` and ``layout.sizes`` do not exist. Mutations (described for the reviewers): make ``offset`` an init field again (the
explicit-offset arm reddens); start the second block one early (the tiling arm reddens). Exact on any host: integer arithmetic.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = [pytest.mark.foundation, pytest.mark.rests_on("tests/gates/numerics/w1/test_face_measure_home.py")]

#: Ordered ``(key, shape)`` sequences; sizes include 1 (a block of one coordinate is where an off-by-one start hides).
SEQUENCES = (
    (("a", (1,)),),
    (("a", (1,)), ("b", (3, 2)), ("c", (1, 1)), ("d", (4,))),
    (("x", (2, 3)), ("y", (1,)), ("z", (5,))),
    ((("lvl", 1, "p"), (2,)), (("lvl", 1, "q"), (1,)), (("lvl", 2, "p"), (3, 1))),
)


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def hand_offsets(pairs: tuple) -> tuple[int, ...]:
    sizes = [int(np.prod(shape)) for _, shape in pairs]
    return tuple(int(sum(sizes[:i])) for i in range(len(sizes)))


def test_offsets_tile_exactly() -> None:
    from orpheus.numerics.coproduct_layout import CoproductLayout, LayoutBlock

    for pairs in SEQUENCES:
        layout = CoproductLayout.from_named_shapes(pairs)
        n = sum(int(np.prod(shape)) for _, shape in pairs)
        _require(layout.total_size == n, f"{pairs}: total_size {layout.total_size}, expected {n}")
        blocks = list(layout.blocks.values())
        _require([b.key for b in blocks] == [k for k, _ in pairs], f"{pairs}: the blocks are not in the order of the pairs")
        _require(tuple(b.offset for b in blocks) == hand_offsets(pairs), f"{pairs}: offsets {[b.offset for b in blocks]} are not the cumulative sizes {hand_offsets(pairs)}")
        _require(tuple(layout.offsets) == hand_offsets(pairs), f"{pairs}: layout.offsets {tuple(layout.offsets)} is not {hand_offsets(pairs)}")
        _require(tuple(layout.sizes) == tuple(int(np.prod(s)) for _, s in pairs), f"{pairs}: layout.sizes {tuple(layout.sizes)}")
        cover = np.zeros(n, dtype=int)
        for b in blocks:
            cover[b.offset : b.offset + b.flat_size] += 1
        _require(bool(np.all(cover == 1)), f"{pairs}: the blocks do not tile range({n}) exactly (coverage counts {np.bincount(cover).tolist()})")
    with pytest.raises(TypeError, match="offset"):
        LayoutBlock(key="a", offset=0, shape=(2,))  # type: ignore[call-arg]
    with pytest.raises(TypeError, match="flat_size"):
        LayoutBlock(key="a", shape=(2,), flat_size=2)  # type: ignore[call-arg]
    _require(LayoutBlock(key="a", shape=(2, 3)).flat_size == 6, "flat_size is not the product of the shape")
