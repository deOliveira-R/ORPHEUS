r"""Codim-1 face primitives: the flat-buffer layout.

L1 primitive (mathematics, knows no neutrons) shared by every codim-1
(face / edge) function space: :class:`LayoutBlock` / :class:`CoproductLayout` — how a
structured collection of faces is laid out in a single flat backing buffer:
each named face is mapped to a contiguous slice + per-face reshape. The
streaming-normal measure lives in :mod:`orpheus.numerics.face_measure`.

The face NAMES these slots are keyed by (the axis letters and the bijection
``(axis, outward sign) <-> "{axis}{min|max}"``) are geometry and live in
:mod:`orpheus.geometry.faces`; until layer DAG WP-B (2026-10-22) they were the
first half of this module. Neither concept here reads them: a layout is
key-generic, and the measure takes the normal as an argument.

Motivation (Depth B step D-G)
==============================

Pre-D-G ``orpheus.sn.boundary_flux.AngularBoundaryFlux`` carried a
dict-like collection of per-geometry-conditional ndarray attributes
(``xmin_face`` / ``xmax_face`` for 1-D, ``xmin_xmax_buf`` /
``ymin_ymax_buf`` for 2-D). Arithmetic dispatched per-attribute with
None-check branches. Cross-method generality (MoC's per-track-family
"faces" in the thousands) would have made this dispatch quadratic in
Python overhead.

Post-D-G :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
inherits :class:`~orpheus.numerics.field.Field` and stores values in a
SINGLE flat backing buffer. :class:`CoproductLayout` is the descriptor that
recovers per-face access via O(1) slice views (memory-shared with the
backing buffer; no copies). Field-inherited algebra operates on the
flat buffer in ONE numpy call regardless of face count — O(1) Python
overhead per arithmetic operation.

The descriptor is geometry-agnostic. The transport angular trace's
:func:`~orpheus.transport.phase_spaces.angular_face_layout` (which SN's
:attr:`~orpheus.sn.problem.SNProblem.angular_trace` carries as its ``layout``) is the
Cartesian-face layout (1-4 faces). A future MoC mesh would return a
per-track-family layout (thousands of faces). The same flat-buffer
arithmetic works for both.

References
==========

* Depth B plan §3.4 (Option Ω flat-buffer storage).
* Plan §6 step D-G (AngularBoundaryFlux as pure Field + SweepScratch carve).
* ``coding-elegance`` Pattern 5 (build the right primitive, not the
  right product) — CoproductLayout is the primitive; the per-method
  AngularBoundaryFlux constructions are the products.
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
    flat_size : int
        Size of the block's flat region. Equals ``prod(shape)``.

    Attributes
    ----------
    offset : int
        Start index of the block's flat region in the backing buffer. It is
        ``init=False``: the owning layout derives it as the sum of the sizes
        of the blocks before this one, so a block whose offset disagrees with
        its position is unspellable (law L13: offsets are computed once).

    Notes
    -----
    The invariant ``flat_size == prod(shape)`` is enforced at construction.
    """

    key: K
    shape: tuple[int, ...]
    flat_size: int
    offset: int = field(init=False)

    def __post_init__(self) -> None:
        expected = int(np.prod(self.shape))
        if self.flat_size != expected:
            raise ValueError(
                f"LayoutBlock({self.key!r}): flat_size={self.flat_size} does "
                f"not match prod(shape={self.shape!r}) = {expected}"
            )

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
            size = int(np.prod(shape))
            block = LayoutBlock(key=key, shape=shape, flat_size=size)
            object.__setattr__(block, "offset", cursor)
            blocks[key] = block
            cursor += size
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

    def pack(self, face_arrays: Mapping[K, NDArray]) -> NDArray:
        r"""Pack per-block ndarrays into a fresh flat backing buffer.

        ``face_arrays`` must cover EVERY block of this layout, each array
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
        provided = set(face_arrays.keys())
        expected = set(self.blocks.keys())
        if provided != expected:
            raise ValueError(
                f"CoproductLayout.pack: face_arrays keys "
                f"{sorted(provided, key=repr)!r} do not match the "
                f"layout's blocks {sorted(expected, key=repr)!r}"
            )
        arrays = {key: np.asarray(face_arrays[key]) for key in self.blocks}
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

    @classmethod
    def from_named_shapes(
        cls, named_shapes: "Sequence[tuple[K, tuple[int, ...]]]",
    ) -> "CoproductLayout[K]":
        r"""Build a :class:`CoproductLayout` from an ordered list of ``(key, shape)`` pairs.

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
