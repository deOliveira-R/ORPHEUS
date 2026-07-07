r"""Codim-1 face primitives: flat-buffer layout + the streaming-normal measure.

L1 primitives (mathematics, knows no neutrons) shared by every codim-1
(face / edge) function space. Two concerns live here:

* :class:`FaceSlot` / :class:`FaceLayout` — how a structured collection of
  faces is laid out in a single flat backing buffer: each named face is
  mapped to a contiguous slice + per-face reshape.
* :func:`face_streaming_normal` — the spatial-trace partial-current
  measure :math:`|\Omega\cdot\hat n|\,w` (the magnitude of the normal
  streaming flux through a boundary face), the metric under which the
  boundary operators are self-adjoint. This is the SPATIAL-trace measure
  ONLY; the angular-pole ψ½ block does NOT route through it (its metric is
  the radial cell volume :math:`V_{\rm cell}`, a state metric — ERR-067;
  see :func:`face_streaming_normal`).

Motivation (Depth B step D-G)
==============================

Pre-D-G :class:`~orpheus.sn.boundary_flux.AngularBoundaryFlux` carried a
dict-like collection of per-geometry-conditional ndarray attributes
(``xmin_face`` / ``xmax_face`` for 1-D, ``xmin_xmax_buf`` /
``ymin_ymax_buf`` for 2-D). Arithmetic dispatched per-attribute with
None-check branches. Cross-method generality (MoC's per-track-family
"faces" in the thousands) would have made this dispatch quadratic in
Python overhead.

Post-D-G :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
inherits :class:`~orpheus.numerics.field.Field` and stores values in a
SINGLE flat backing buffer. :class:`FaceLayout` is the descriptor that
recovers per-face access via O(1) slice views (memory-shared with the
backing buffer; no copies). Field-inherited algebra operates on the
flat buffer in ONE numpy call regardless of face count — O(1) Python
overhead per arithmetic operation.

The descriptor is geometry-agnostic. SN's
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout` returns the
Cartesian-face layout (1-4 faces). A future MoC mesh would return a
per-track-family layout (thousands of faces). The same flat-buffer
arithmetic works for both.

References
==========

* Depth B plan §3.4 (Option Ω flat-buffer storage).
* Plan §6 step D-G (AngularBoundaryFlux as pure Field + SweepScratch carve).
* ``coding-elegance`` Pattern 5 (build the right primitive, not the
  right product) — FaceLayout is the primitive; the per-method
  AngularBoundaryFlux constructions are the products.
"""

from __future__ import annotations

from collections.abc import Hashable
from dataclasses import dataclass
from typing import Generic, Mapping, Sequence, TypeVar

import numpy as np
from numpy.typing import NDArray


__all__ = ["AXIS_NAMES", "FaceSlot", "FaceLayout", "face_streaming_normal"]

#: The slot-key type of a :class:`FaceLayout` / :class:`FaceSlot`. The
#: descriptor is key-generic: SN Cartesian spatial traces key faces by
#: ``str`` name (``FaceLayout[str]``); the ψ½ starting-direction space keys
#: each ``(level, sign)`` leg's ``cells`` / ``corner`` slots by the
#: structured tuple ``(level, sign, part)``
#: (``FaceLayout[tuple[int, int, str]]``). The flat-buffer discipline
#: (offset assignment, slice views, gap/overlap validation) is identical
#: across key regimes — only the key's *type* differs. Bounded by
#: :class:`~collections.abc.Hashable`: the key is a mapping key, so an
#: unhashable realization (``FaceLayout[list[int]]``) is unspellable.
K = TypeVar("K", bound=Hashable)


#: Spatial axis names, positional-by-axis — the single source of the
#: axis↔name crosswalk for the ``"{axis}min"`` / ``"{axis}max"``
#: boundary-face naming convention that :class:`FaceLayout`, the trace
#: space, :attr:`SNMesh.bc <orpheus.sn.mesh.augmented_mesh.SNMesh.bc>`, and the
#: sweep schedule all key on. No consumer hand-lists ``("x", ...),
#: ("y", ...)`` pairs — every face-name derivation
#: (:attr:`FaceLabel.face_name <orpheus.transport.mesh.axis.FaceLabel.face_name>`,
#: the walk's in/outflow faces, the schedule's outgoing faces, the
#: trace's outward-normal table) renders through this tuple. Lives
#: here, at the bottom of the dependency graph next to
#: :class:`FaceLayout` (the keeper of the face-name string world);
#: :mod:`orpheus.transport.mesh.axis` imports it upward (moved from there in C5.3,
#: #225, so the geometry-blind trace space could share it without an
#: sn-ward import).
AXIS_NAMES = ("x", "y", "z")


@dataclass(frozen=True)
class FaceSlot(Generic[K]):
    r"""Per-face slot in a flat boundary buffer.

    Parameters
    ----------
    key : K
        The slot's identity in the owning :class:`FaceLayout`.
        A ``str`` face name (``"xmin"`` / ``"xmax"`` / ``"ymin"`` / ``"ymax"``
        for SN Cartesian spatial traces; arbitrary tokens for other methods)
        in the ``FaceLayout[str]`` case, or a structured key when the faces
        are not stringly named — the ψ½ starting-direction space keys each
        slot by ``(level, sign, part)`` (a ``FaceLayout[tuple[int, int, str]]``).
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
                f"FaceSlot({self.key!r}): flat_size={self.flat_size} does "
                f"not match prod(shape={self.shape!r}) = {expected}"
            )
        if self.offset < 0:
            raise ValueError(
                f"FaceSlot({self.key!r}): offset must be non-negative; "
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
class FaceLayout(Generic[K]):
    r"""Descriptor mapping face keys to slices in a flat backing buffer.

    Parameters
    ----------
    faces : Mapping[K, FaceSlot[K]]
        Ordered mapping from face key to :class:`FaceSlot`. Iteration
        order determines flat-buffer concatenation order (which matters
        for reproducible serialization and snapshot comparison).
    total_size : int
        Total length of the flat backing buffer. Must equal the sum of
        each face's :attr:`FaceSlot.flat_size`.

    Notes
    -----
    The descriptor is IMMUTABLE; the backing buffer it describes is
    owned elsewhere. A single :class:`FaceLayout` can be shared across
    arbitrarily many :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
    instances on the same mesh (validity invariant: same layout
    identity ⟹ compatible algebra).
    """

    faces: Mapping[K, FaceSlot[K]]
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
                    f"FaceLayout: slot {slot.key!r} starts at offset "
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
        cls, named_shapes: "Sequence[tuple[K, tuple[int, ...]]]",
    ) -> "FaceLayout[K]":
        r"""Build a :class:`FaceLayout` from an ordered list of ``(key, shape)`` pairs.

        Iteration order of ``named_shapes`` determines the flat-buffer
        offset assignment (first face gets offset 0, second face gets
        offset = first face's flat_size, etc.). Use this constructor
        when the per-method mesh has a natural face ordering it wants
        to expose. Key-generic: the SN spatial trace passes ``str`` face
        names (``FaceLayout[str]``); the ψ½ starting-direction space passes
        ``(level, sign, part)`` tuples (``FaceLayout[tuple[int, int, str]]``).

        Parameters
        ----------
        named_shapes : Sequence[tuple[K, tuple[int, ...]]]
            Ordered list of ``(face_key, per_face_shape)`` pairs.

        Returns
        -------
        FaceLayout[K]
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
        faces: dict[K, FaceSlot[K]] = {}
        cursor = 0
        for key, shape in named_shapes:
            flat_size = int(np.prod(shape))
            faces[key] = FaceSlot(
                key=key, offset=cursor, shape=shape, flat_size=flat_size,
            )
            cursor += flat_size
        return cls(faces=faces, total_size=cursor)


# ─────────────────────────────────────────────────────────────────────
# The codim-1 face measure
# ─────────────────────────────────────────────────────────────────────


def face_streaming_normal(
    normal_coefficient: NDArray, quadrature_weight: NDArray,
) -> NDArray:
    r"""The spatial-trace partial-current measure :math:`|\Omega\cdot\hat n_f|\,w`.

    The magnitude of the *normal component of the streaming flux through a
    codim-1 spatial boundary face*, weighted by the angular quadrature:
    with :math:`c = \Omega\cdot\hat n_f` it is the partial-current metric
    :math:`G_s = |\Omega\cdot\hat n_f|\,w_n` (Lewis & Miller §3.7), the
    boundary inner product under which the
    :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
    boundary operators are self-adjoint. It vanishes exactly at
    **grazing** incidence (:math:`\Omega \perp \hat n_f`, so
    :math:`|\Omega\cdot\hat n_f| = 0`) — a tangential ordinate carries no
    partial current across the face.

    L1 primitive (mathematics, knows no neutrons). Its single production
    consumer is the trace-metric build in
    :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`.

    .. note::

       **This is the SPATIAL-trace measure only — it is NOT the ψ½
       starting-direction (angular-pole) block's metric.**  An earlier
       design (``facefield_codim1_design.md`` §3.2) proposed unifying both
       codim-1 metrics through this one kernel, reading the pole's angular
       through-flux coefficient :math:`(1-\mu^2)\big|_{\mu=\pm1} = 0` as
       the ψ½ block's Hilbert inner product.  That unification is
       **refuted** (ERR-067): :math:`(1-\mu^2)|_{\rm pole}` is an
       *operator coefficient* (the angular-redistribution strength that
       makes :math:`\mu=\pm1` a straight characteristic), **not** a state
       metric.  ψ½ is a first-class radial state field whose Hilbert
       metric is the radial cell volume :math:`G_{\rm sd} = V_{\rm cell}`
       (SPD, nonzero) — set by its operator role, not by an integration
       weight.  The two codim-1 metrics do **not** share a kernel; the
       metric is a per-leaf property (spatial trace → this partial
       current; pole → :math:`V_{\rm cell}`), exactly as the bulk's metric
       is its own per-leaf :math:`V\,w`.  See
       :mod:`orpheus.numerics.spaces.radial_characteristic_space`.

    Parameters
    ----------
    normal_coefficient : NDArray
        The ordinate's signed normal streaming coefficient
        :math:`c = \Omega\cdot\hat n_f` on a spatial boundary face. Any
        shape.
    quadrature_weight : NDArray
        The angular quadrature weight :math:`w`, broadcast-compatible with
        ``normal_coefficient``.

    Returns
    -------
    NDArray
        :math:`|c|\,w`, the elementwise normal streaming measure.

    References
    ----------
    * Lewis, E.E. & Miller, W.F. (1993). *Computational Methods of Neutron
      Transport*. §3.7 (the partial-current boundary metric).
    """
    return np.abs(normal_coefficient) * quadrature_weight
