r"""The coproduct layout: how the blocks of a direct sum, or the faces of a boundary trace, sit in one flat backing buffer.

L1 primitive (mathematics, knows no neutrons). :class:`CoproductLayout` partitions one flat buffer into keyed blocks in order: each key owns a
contiguous region of the leading axis and a shape that region is reshaped to. The coproduct reading is that the buffer is the direct sum of the blocks,
:meth:`~CoproductLayout.pack` is the sum of the injections, :meth:`~CoproductLayout.unpack` the projections (views, no copies), and the offsets are
derived once, here, from the sizes (law L13), so no consumer slices a block by hand: it reads :attr:`LayoutBlock.region`,
:meth:`~CoproductLayout.slice_of` or :meth:`LayoutBlock.slice_view`.

Two families of consumers share it. The summands of a :class:`~orpheus.numerics.direct_sum.DirectSumSpace` are keyed by their labels. The faces
of a boundary trace space are keyed by name (the axis letters and the bijection ``(axis, outward sign) <-> "{axis}{min|max}"`` are geometry and live in
:mod:`orpheus.geometry.faces`; a layout is key-generic and reads none of it) or, for the ψ½ starting-direction space, by ``(level, sign, part)``.
The streaming-normal measure lives in :mod:`orpheus.numerics.face_measure`.

Why one flat buffer. A boundary field stores its values in a single flat backing buffer
(:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` inherits :class:`~orpheus.numerics.field.Field`), so the field algebra
is one numpy call whatever the number of blocks, where a per-face attribute dispatch grew quadratic with the face count. The layout is the descriptor that
recovers per-block access as views. It is geometry-agnostic: the Cartesian layout has one to four faces and a method-of-characteristics mesh would
have thousands.

References
==========

* Depth B plan 3.4 (option Omega, flat-buffer storage) and plan 6 step D-G (the boundary flux as a pure ``Field``).
* ``coding-elegance`` Pattern 5 (build the right primitive, not the right product): the layout is the primitive and the per-method boundary fields are
  the products.
"""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass, field
from typing import Generic, Mapping, Sequence, TypeVar

import numpy as np
from numpy.typing import NDArray


__all__ = [
    "LayoutBlock",
    "CoproductLayout",
]

#: The slot-key type of a :class:`CoproductLayout` / :class:`LayoutBlock`. The
#: descriptor is key-generic: SN Cartesian spatial traces key faces by
#: ``str`` name (``CoproductLayout[str]``); the ψ½ starting-direction space keys
#: each ``(level, sign)`` leg's ``cells`` / ``corner`` slots by the
#: structured tuple ``(level, sign, part)``
#: (``CoproductLayout[tuple[int, int, str]]``). The flat-buffer discipline
#: (offset derivation, slice views) is identical
#: across key regimes — only the key's *type* differs. Bounded by
#: :class:`~collections.abc.Hashable`: the key is a mapping key, so an
#: unhashable realization (``CoproductLayout[list[int]]``) is unspellable.
K = TypeVar("K", bound=Hashable)


@dataclass(frozen=True)
class LayoutBlock(Generic[K]):  # concept: layout | role: part | kind: set
    r"""One block of a :class:`CoproductLayout`: a key, a shape and the flat region it owns.

    Parameters
    ----------
    key : K
        The block's identity in the owning :class:`CoproductLayout`: a ``str``
        face name (``"xmin"`` / ``"xmax"`` / ``"ymin"`` / ``"ymax"`` for SN
        Cartesian spatial traces), a summand label of a direct sum, or a
        structured key when the blocks are not stringly named (the ψ½
        starting-direction space keys each slot by ``(level, sign, part)``).
    shape : tuple[int, ...]
        The block's shape. The flat region is reshaped to this shape (plus the
        trailing payload axes of the buffer) when sliced.

    Attributes
    ----------
    flat_size : int
        Size of the block's flat region, ``prod(shape)``: a property over the shape, so no copy of it can disagree.
    offset : int
        Start index of the block's flat region in the backing buffer. It is
        ``init=False``: the owning layout derives it as the sum of the sizes
        of the blocks before this one, so a block whose offset disagrees with
        its position is unspellable (law L13: offsets are computed once).

    """

    key: K
    shape: tuple[int, ...]
    offset: int = field(init=False)

    @property
    def flat_size(self) -> int:
        r"""The size of the block's flat region: the product of the shape."""
        return int(np.prod(self.shape))

    @property
    def region(self) -> slice:
        r"""The block's flat region of the leading axis."""
        return slice(self.offset, self.offset + self.flat_size)

    def slice_view(self, flat_buf: NDArray) -> NDArray:
        r"""Return a shaped view into the flat buffer — no copy.

        The buffer's leading axis is the layout's; its trailing axes (a payload,
        ``flat_buf.shape[1:]``) are kept after the block's own shape. The
        returned array shares memory with ``flat_buf``; writes to the view
        propagate to the backing buffer (the intended access pattern for
        sweep-side face writes that must be visible to the BC partner on the
        next iteration).
        """
        return flat_buf[self.region].reshape(self.shape + tuple(flat_buf.shape[1:]))


@dataclass(frozen=True)
class CoproductLayout(Generic[K]):  # concept: layout | role: arm | kind: set
    r"""The partition of one flat backing buffer into keyed blocks, in order.

    Parameters
    ----------
    named_shapes : tuple[tuple[K, tuple[int, ...]], ...]
        Ordered ``(key, shape)`` pairs. The order is the flat concatenation
        order (which matters for reproducible serialization and snapshot
        comparison) and the keys are pairwise distinct.

    Attributes
    ----------
    blocks : Mapping[K, LayoutBlock[K]]
        The blocks by key, in buffer order; each block's offset is the sum of
        the sizes before it (derived, so the blocks tile ``range(total_size)``
        exactly).
    total_size : int
        Length of the flat backing buffer: the sum of the block sizes.

    Notes
    -----
    The descriptor is IMMUTABLE; the backing buffer it describes is owned
    elsewhere. A single :class:`CoproductLayout` can be shared across
    arbitrarily many fields on the same mesh. It is the one place the offsets
    of a direct sum are computed (law L13): the summands of a
    :class:`~orpheus.numerics.direct_sum.DirectSumSpace` read their slices
    from here.
    """

    named_shapes: tuple[tuple[K, tuple[int, ...]], ...]
    blocks: Mapping[K, LayoutBlock[K]] = field(init=False, compare=False, repr=False)
    total_size: int = field(init=False)

    def __post_init__(self) -> None:
        pairs = tuple((key, tuple(int(n) for n in shape)) for key, shape in self.named_shapes)
        object.__setattr__(self, "named_shapes", pairs)
        keys = [key for key, _ in pairs]
        if len(set(keys)) != len(keys):
            raise ValueError(
                f"CoproductLayout: block keys must be distinct; got {sorted(keys, key=repr)!r}"
            )
        blocks: dict[K, LayoutBlock[K]] = {}
        cursor = 0
        for key, shape in pairs:
            block = LayoutBlock(key=key, shape=shape)
            object.__setattr__(block, "offset", cursor)
            blocks[key] = block
            cursor += block.flat_size
        object.__setattr__(self, "blocks", blocks)
        object.__setattr__(self, "total_size", cursor)

    @property
    def sizes(self) -> tuple[int, ...]:
        r"""The flat size of each block, in buffer order."""
        return tuple(block.flat_size for block in self.blocks.values())

    @property
    def offsets(self) -> tuple[int, ...]:
        r"""The start of each block's flat region, in buffer order."""
        return tuple(block.offset for block in self.blocks.values())

    def slice_of(self, key: K) -> slice:
        r"""The slice of the flat buffer's leading axis that holds block ``key``.

        Raises
        ------
        KeyError
            If ``key`` is not a block of this layout (the message lists the keys).
        """
        if key not in self.blocks:
            raise KeyError(
                f"CoproductLayout.slice_of: {key!r} is not a block; blocks are {list(self.blocks)!r}"
            )
        return self.blocks[key].region

    @property
    def structure(self) -> tuple[tuple[K, int, int], ...]:
        r"""``(key, offset, flat size)`` per slot, in buffer order.

        The layout's structural content: what a trace space's name digest
        reads through :func:`~orpheus.values.content.name_digest` (#405 P1
        step 5). A tuple, because the slot ORDER is the buffer's order and so
        is content, where a mapping's order would not be.
        """
        return tuple(
            (key, int(block.offset), int(block.flat_size))
            for key, block in self.blocks.items()
        )

    def pack(self, parts: Mapping[K, NDArray]) -> NDArray:
        r"""Pack per-block ndarrays into a fresh flat backing buffer: the injections summed.

        ``parts`` must cover EVERY block of this layout, each array
        matching its block's shape followed by one common trailing payload
        (none for a plain buffer; ``G (x) I_P`` on a payload of ``P`` columns).
        The packing loop lives HERE because every fact it reads — the block
        set, each block's shape, the slice map, the total size — is layout
        knowledge (native place, CS4b S6.2); the field factories
        (:meth:`~orpheus.transport.fields._bases.BoundaryField.from_face_arrays`)
        are typed entries over this.

        Raises
        ------
        ValueError
            If the keys differ from this layout's blocks, any array's shape
            mismatches its block, or the payloads of two blocks differ.
        """
        provided = set(parts.keys())
        expected = set(self.blocks.keys())
        if provided != expected:
            raise ValueError(
                f"CoproductLayout.pack: parts keys "
                f"{sorted(provided, key=repr)!r} do not match the "
                f"layout's blocks {sorted(expected, key=repr)!r}"
            )
        arrays = {key: np.asarray(parts[key]) for key in self.blocks}
        lead = next(iter(self.blocks.items()), None)
        payload = () if lead is None else tuple(arrays[lead[0]].shape[len(lead[1].shape):])
        flat = np.zeros((self.total_size,) + payload)
        for name, block in self.blocks.items():
            arr = arrays[name]
            if arr.shape != block.shape + payload:
                raise ValueError(
                    f"CoproductLayout.pack: block {name!r} array shape "
                    f"{arr.shape!r} does not match its block shape "
                    f"{block.shape!r} with the common payload {payload!r}"
                )
            block.slice_view(flat)[:] = arr
        return flat

    def unpack(self, flat: NDArray) -> dict[K, NDArray]:
        r"""The per-block views of a flat buffer: the inverse of :meth:`pack`.

        ``flat`` has this layout's ``total_size`` on its leading axis and any
        trailing payload; each block comes back with its own shape followed by
        that payload, as a view that shares memory with ``flat``.

        Raises
        ------
        ValueError
            If ``flat`` is a scalar or its leading axis is not ``total_size``.
        """
        buf = np.asarray(flat)
        if buf.ndim < 1 or buf.shape[0] != self.total_size:
            raise ValueError(
                f"CoproductLayout.unpack: the leading axis of the buffer must have "
                f"total_size={self.total_size}; got shape {buf.shape!r}"
            )
        return {key: block.slice_view(buf) for key, block in self.blocks.items()}

    def axis_index(self, axis: int) -> NDArray:
        r"""The index along block axis ``axis`` of every flat coordinate, in the flat chart's order (an integer array of ``total_size``).

        For each block the index runs along that axis of the block's shape and is constant along the others, so a block ``(2, ng, n)`` with
        ``axis=1`` gives the group of each of its coordinates. Every block must have that axis. The layout computes this because the flat region of
        each block is its own knowledge (law L13); a consumer reads it and never slices a block by hand.
        """
        out = np.zeros(self.total_size, dtype=np.intp)
        for block in self.blocks.values():
            if not -len(block.shape) <= axis < len(block.shape):
                raise ValueError(f"CoproductLayout.axis_index: block {block.key!r} of shape {block.shape} has no axis {axis}")
            ax = axis % len(block.shape)
            index = np.arange(block.shape[ax]).reshape(tuple(-1 if i == ax else 1 for i in range(len(block.shape))))
            out[block.region] = np.broadcast_to(index, block.shape).ravel()
        return out

    @classmethod
    def from_named_shapes(
        cls, named_shapes: "Sequence[tuple[K, tuple[int, ...]]]",
    ) -> "CoproductLayout[K]":
        r"""Build a :class:`CoproductLayout` from an ordered list of ``(key, shape)`` pairs.

        The constructor takes the same pairs as ``named_shapes=``; this door stays because its 56 callers read it as a named construction of a layout
        from an ordered list, and a second spelling of the constructor's one argument carries no second behaviour.

        Iteration order of ``named_shapes`` determines the flat-buffer
        offset assignment (first block gets offset 0, second block gets
        offset = first block's flat_size, etc.). Use this constructor
        when the per-method mesh has a natural face ordering it wants
        to expose. Key-generic: the SN spatial trace passes ``str`` face
        names (``CoproductLayout[str]``); the ψ½ starting-direction space passes
        ``(level, sign, part)`` tuples (``CoproductLayout[tuple[int, int, str]]``).

        Examples
        --------
        >>> layout = CoproductLayout.from_named_shapes([
        ...     ("xmin", (2, 1)),
        ...     ("xmax", (2, 1)),
        ... ])
        >>> layout.total_size
        4
        >>> layout.blocks["xmin"].offset
        0
        >>> layout.blocks["xmax"].offset
        2
        """
        return cls(named_shapes=tuple((key, tuple(shape)) for key, shape in named_shapes))
