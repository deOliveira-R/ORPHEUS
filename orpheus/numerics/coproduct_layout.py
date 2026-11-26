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
from dataclasses import dataclass
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
#: (offset assignment, slice views, gap/overlap validation) is identical
#: across key regimes — only the key's *type* differs. Bounded by
#: :class:`~collections.abc.Hashable`: the key is a mapping key, so an
#: unhashable realization (``CoproductLayout[list[int]]``) is unspellable.
K = TypeVar("K", bound=Hashable)


@dataclass(frozen=True)
class LayoutBlock(Generic[K]):  # concept: layout | role: part | kind: set
    r"""Per-face slot in a flat boundary buffer.

    Parameters
    ----------
    key : K
        The slot's identity in the owning :class:`CoproductLayout`.
        A ``str`` face name (``"xmin"`` / ``"xmax"`` / ``"ymin"`` / ``"ymax"``
        for SN Cartesian spatial traces; arbitrary tokens for other methods)
        in the ``CoproductLayout[str]`` case, or a structured key when the faces
        are not stringly named — the ψ½ starting-direction space keys each
        slot by ``(level, sign, part)`` (a ``CoproductLayout[tuple[int, int, str]]``).
    offset : int
        Start index of this face's flat region in the backing buffer.
    shape : tuple[int, ...]
        Per-face shape. The flat region is reshaped to this shape when
        sliced (e.g. ``(N, ng)`` for a 1-D slab face, ``(N, ng, ny)`` for
        a 2-D Cartesian x-face).
    flat_size : int
        Size of this face's flat region. Must equal ``prod(shape)``.

    Notes
    -----
    The validity invariant ``flat_size == prod(shape)`` is enforced at
    construction. Consumers may rely on this without re-checking.
    """

    key: K
    offset: int
    shape: tuple[int, ...]
    flat_size: int

    def __post_init__(self) -> None:
        expected = int(np.prod(self.shape))
        if self.flat_size != expected:
            raise ValueError(
                f"LayoutBlock({self.key!r}): flat_size={self.flat_size} does "
                f"not match prod(shape={self.shape!r}) = {expected}"
            )
        if self.offset < 0:
            raise ValueError(
                f"LayoutBlock({self.key!r}): offset must be non-negative; "
                f"got {self.offset}"
            )

    def slice_view(self, flat_buf: NDArray) -> NDArray:
        r"""Return a shaped view into the flat buffer — no copy.

        The returned array shares memory with ``flat_buf``; writes to
        the view propagate to the backing buffer (the intended access
        pattern for sweep-side face writes that must be visible to the
        BC partner on the next iteration).
        """
        return flat_buf[self.offset : self.offset + self.flat_size].reshape(self.shape)


@dataclass(frozen=True)
class CoproductLayout(Generic[K]):  # concept: layout | role: arm | kind: set
    r"""Descriptor mapping face keys to slices in a flat backing buffer.

    Parameters
    ----------
    blocks : Mapping[K, LayoutBlock[K]]
        Ordered mapping from face key to :class:`LayoutBlock`. Iteration
        order determines flat-buffer concatenation order (which matters
        for reproducible serialization and snapshot comparison).
    total_size : int
        Total length of the flat backing buffer. Must equal the sum of
        each face's :attr:`LayoutBlock.flat_size`.

    Notes
    -----
    The descriptor is IMMUTABLE; the backing buffer it describes is
    owned elsewhere. A single :class:`CoproductLayout` can be shared across
    arbitrarily many :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
    instances on the same mesh (validity invariant: same layout
    identity ⟹ compatible algebra).
    """

    blocks: Mapping[K, LayoutBlock[K]]
    total_size: int

    def __post_init__(self) -> None:
        expected = sum(slot.flat_size for slot in self.blocks.values())
        if self.total_size != expected:
            raise ValueError(
                f"CoproductLayout: total_size={self.total_size} does not match "
                f"sum(slot.flat_size for face slots) = {expected}"
            )
        # Validate no slot overlap or gap.
        sorted_slots = sorted(self.blocks.values(), key=lambda s: s.offset)
        cursor = 0
        for slot in sorted_slots:
            if slot.offset != cursor:
                raise ValueError(
                    f"CoproductLayout: slot {slot.key!r} starts at offset "
                    f"{slot.offset} but cursor is at {cursor} "
                    f"(gap or overlap with previous slot)"
                )
            cursor += slot.flat_size
        if cursor != self.total_size:
            raise ValueError(
                f"CoproductLayout: cursor reached {cursor} after walking all "
                f"slots but total_size={self.total_size}"
            )

    @property
    def structure(self) -> tuple[tuple[K, int, int], ...]:
        r"""``(key, offset, flat size)`` per slot, in buffer order.

        The layout's structural content: what a trace space's name digest
        reads through :func:`~orpheus.values.content.name_digest` (#405 P1
        step 5). A tuple, because the slot ORDER is the buffer's order and so
        is content, where a mapping's order would not be.
        """
        return tuple(
            (key, int(slot.offset), int(slot.flat_size))
            for key, slot in self.blocks.items()
        )

    def pack(self, face_arrays: Mapping[K, NDArray]) -> NDArray:
        r"""Pack per-face ndarrays into a fresh flat backing buffer.

        ``face_arrays`` must cover EVERY face of this layout, each array
        matching its slot's shape. The packing loop lives HERE because
        every fact it reads — the face set, each slot's shape, the slice
        map, the total size — is layout knowledge (native place, CS4b
        S6.2); the field factories
        (:meth:`~orpheus.transport.fields._bases.BoundaryField.from_face_arrays`)
        are typed entries over this.

        Raises
        ------
        ValueError
            If the keys differ from this layout's faces, or any per-face
            array's shape mismatches its slot.
        """
        provided = set(face_arrays.keys())
        expected = set(self.blocks.keys())
        if provided != expected:
            raise ValueError(
                f"CoproductLayout.pack: face_arrays keys "
                f"{sorted(provided, key=repr)!r} do not match the "
                f"layout's faces {sorted(expected, key=repr)!r}"
            )
        flat = np.zeros(self.total_size)
        for name, slot in self.blocks.items():
            arr = np.asarray(face_arrays[name])
            if arr.shape != slot.shape:
                raise ValueError(
                    f"CoproductLayout.pack: face {name!r} array shape "
                    f"{arr.shape!r} does not match its slot shape "
                    f"{slot.shape!r}"
                )
            slot.slice_view(flat)[:] = arr
        return flat

    @classmethod
    def from_named_shapes(
        cls, named_shapes: "Sequence[tuple[K, tuple[int, ...]]]",
    ) -> "CoproductLayout[K]":
        r"""Build a :class:`CoproductLayout` from an ordered list of ``(key, shape)`` pairs.

        Iteration order of ``named_shapes`` determines the flat-buffer
        offset assignment (first face gets offset 0, second face gets
        offset = first face's flat_size, etc.). Use this constructor
        when the per-method mesh has a natural face ordering it wants
        to expose. Key-generic: the SN spatial trace passes ``str`` face
        names (``CoproductLayout[str]``); the ψ½ starting-direction space passes
        ``(level, sign, part)`` tuples (``CoproductLayout[tuple[int, int, str]]``).

        Parameters
        ----------
        named_shapes : Sequence[tuple[K, tuple[int, ...]]]
            Ordered list of ``(face_key, per_face_shape)`` pairs.

        Returns
        -------
        CoproductLayout[K]
            The constructed layout.

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
        faces: dict[K, LayoutBlock[K]] = {}
        cursor = 0
        for key, shape in named_shapes:
            flat_size = int(np.prod(shape))
            faces[key] = LayoutBlock(
                key=key, offset=cursor, shape=shape, flat_size=flat_size,
            )
            cursor += flat_size
        return cls(blocks=faces, total_size=cursor)
