r"""Face layout descriptor for flat boundary-trace buffers.

L1 primitive (mathematics, knows no neutrons). Describes how a structured
collection of boundary faces is laid out in a single flat backing buffer:
each named face is mapped to a contiguous slice + per-face reshape.

Motivation (Depth B step D-G)
==============================

Pre-D-G :class:`~orpheus.sn.boundary_flux.BoundaryFlux` carried a
dict-like collection of per-geometry-conditional ndarray attributes
(``xmin_face`` / ``xmax_face`` for 1-D, ``xmin_xmax_buf`` /
``ymin_ymax_buf`` for 2-D). Arithmetic dispatched per-attribute with
None-check branches. Cross-method generality (MoC's per-track-family
"faces" in the thousands) would have made this dispatch quadratic in
Python overhead.

Post-D-G :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
inherits :class:`~orpheus.numerics.field.Field` and stores values in a
SINGLE flat backing buffer. :class:`FaceLayout` is the descriptor that
recovers per-face access via O(1) slice views (memory-shared with the
backing buffer; no copies). Field-inherited algebra operates on the
flat buffer in ONE numpy call regardless of face count — O(1) Python
overhead per arithmetic operation.

The descriptor is geometry-agnostic. SN's
:meth:`~orpheus.sn.geometry.SNMesh.boundary_face_layout` returns the
Cartesian-face layout (1-4 faces). A future MoC mesh would return a
per-track-family layout (thousands of faces). The same flat-buffer
arithmetic works for both.

References
==========

* Depth B plan §3.4 (Option Ω flat-buffer storage).
* Plan §6 step D-G (BoundaryFlux as pure Field + SweepScratch carve).
* ``coding-elegance`` Pattern 5 (build the right primitive, not the
  right product) — FaceLayout is the primitive; the per-method
  BoundaryFlux constructions are the products.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from numpy.typing import NDArray


__all__ = ["FaceSlot", "FaceLayout"]


@dataclass(frozen=True)
class FaceSlot:
    r"""Per-face slot in a flat boundary buffer.

    Parameters
    ----------
    name : str
        Face name (e.g. ``"xmin"``, ``"xmax"``, ``"ymin"``, ``"ymax"``
        for SN Cartesian geometries; arbitrary tokens for other methods).
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

    name: str
    offset: int
    shape: tuple[int, ...]
    flat_size: int

    def __post_init__(self) -> None:
        expected = int(np.prod(self.shape))
        if self.flat_size != expected:
            raise ValueError(
                f"FaceSlot({self.name!r}): flat_size={self.flat_size} does "
                f"not match prod(shape={self.shape!r}) = {expected}"
            )
        if self.offset < 0:
            raise ValueError(
                f"FaceSlot({self.name!r}): offset must be non-negative; "
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
class FaceLayout:
    r"""Descriptor mapping face names to slices in a flat backing buffer.

    Parameters
    ----------
    faces : Mapping[str, FaceSlot]
        Ordered mapping from face name to :class:`FaceSlot`. Iteration
        order determines flat-buffer concatenation order (which matters
        for reproducible serialization and snapshot comparison).
    total_size : int
        Total length of the flat backing buffer. Must equal the sum of
        each face's :attr:`FaceSlot.flat_size`.

    Notes
    -----
    The descriptor is IMMUTABLE; the backing buffer it describes is
    owned elsewhere. A single :class:`FaceLayout` can be shared across
    arbitrarily many :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
    instances on the same mesh (validity invariant: same layout
    identity ⟹ compatible algebra).
    """

    faces: Mapping[str, FaceSlot]
    total_size: int

    def __post_init__(self) -> None:
        expected = sum(slot.flat_size for slot in self.faces.values())
        if self.total_size != expected:
            raise ValueError(
                f"FaceLayout: total_size={self.total_size} does not match "
                f"sum(slot.flat_size for face slots) = {expected}"
            )
        # Validate no slot overlap or gap.
        sorted_slots = sorted(self.faces.values(), key=lambda s: s.offset)
        cursor = 0
        for slot in sorted_slots:
            if slot.offset != cursor:
                raise ValueError(
                    f"FaceLayout: slot {slot.name!r} starts at offset "
                    f"{slot.offset} but cursor is at {cursor} "
                    f"(gap or overlap with previous slot)"
                )
            cursor += slot.flat_size
        if cursor != self.total_size:
            raise ValueError(
                f"FaceLayout: cursor reached {cursor} after walking all "
                f"slots but total_size={self.total_size}"
            )

    @classmethod
    def from_named_shapes(
        cls, named_shapes: "list[tuple[str, tuple[int, ...]]]",
    ) -> "FaceLayout":
        r"""Build a :class:`FaceLayout` from an ordered list of ``(name, shape)`` pairs.

        Iteration order of ``named_shapes`` determines the flat-buffer
        offset assignment (first face gets offset 0, second face gets
        offset = first face's flat_size, etc.). Use this constructor
        when the per-method mesh has a natural face ordering it wants
        to expose.

        Parameters
        ----------
        named_shapes : list[tuple[str, tuple[int, ...]]]
            Ordered list of ``(face_name, per_face_shape)`` pairs.

        Returns
        -------
        FaceLayout
            The constructed layout.

        Examples
        --------
        >>> layout = FaceLayout.from_named_shapes([
        ...     ("xmin", (2, 1)),
        ...     ("xmax", (2, 1)),
        ... ])
        >>> layout.total_size
        4
        >>> layout.faces["xmin"].offset
        0
        >>> layout.faces["xmax"].offset
        2
        """
        faces: dict[str, FaceSlot] = {}
        cursor = 0
        for name, shape in named_shapes:
            flat_size = int(np.prod(shape))
            faces[name] = FaceSlot(
                name=name, offset=cursor, shape=shape, flat_size=flat_size,
            )
            cursor += flat_size
        return cls(faces=faces, total_size=cursor)
