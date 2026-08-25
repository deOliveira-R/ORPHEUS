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

Pre-D-G ``orpheus.sn.boundary_flux.AngularBoundaryFlux`` carried a
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


__all__ = [
    "AXIS_NAMES",
    "FACE_NAMES",
    "FaceSlot",
    "FaceLayout",
    "face_name",
    "face_normal",
    "face_opposite",
    "face_streaming_normal",
]

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


# ─────────────────────────────────────────────────────────────────────
# The face-name bijection
# ─────────────────────────────────────────────────────────────────────
#
# A boundary face of an axis-aligned mesh is completely described by its
# OUTWARD NORMAL, and for an axis-aligned face that normal is sparse:
# ``n_f = sign · ê_axis``. So ``(axis, outward_sign)`` and the face NAME
# ``"{axis}{min|max}"`` are two encodings of one object, and the map between
# them is a bijection.
#
# It lives here, beside :data:`AXIS_NAMES`, because this module is already the
# keeper of the face-name string world and sits at the bottom of the dependency
# graph (the trace space and :mod:`orpheus.transport.mesh.axis` both import
# upward from it). Before campaign phase **B3.4c** the convention had no single
# home and was transcribed at five sites — ``_FACE_NORMALS`` in the trace space
# (the parse), two verbatim-twin ``f"{axis}{'max' if outward_sign == +1 else
# 'min'}"`` renders in the SN boundary layer, one more in the sweep schedule,
# and a ``face_name.endswith("max")`` reverse parse in the method registry.
# Two of those sites carried a docstring pointing at "the one primitive" that
# did not exist yet. B3.4c needed a SIXTH consumer (the periodic partner face),
# which is what forced the collapse.

#: Outward-normal sign per face-name suffix. The ``min`` face of an axis has
#: the INWARD-pointing unit vector as its outward normal (:math:`-\hat e`),
#: the ``max`` face the outward one — which is the whole content of the
#: ``"{axis}{min|max}"`` convention, and the reason a bare string carries
#: enough orientation to derive :math:`\Omega\cdot\hat n` from.
_FACE_SUFFIX_SIGN: dict[str, int] = {"min": -1, "max": +1}

#: The inverse map, DERIVED rather than transcribed — the two directions of
#: one bijection must not be able to disagree.
_FACE_SIGN_SUFFIX: dict[int, str] = {
    sign: suffix for suffix, sign in _FACE_SUFFIX_SIGN.items()
}


def face_name(axis: int, outward_sign: int) -> str:
    r"""Render ``(axis index, outward-normal sign)`` as the canonical face name.

    The forward direction of the face-name bijection: ``(0, -1) -> "xmin"``.
    Inverse of :func:`face_normal`.

    Parameters
    ----------
    axis : int
        Position of the axis in :data:`AXIS_NAMES`.
    outward_sign : int
        ``+1`` for the axis's ``max`` face, ``-1`` for its ``min`` face.

    Raises
    ------
    ValueError
        If ``axis`` is outside the named-axis inventory, or ``outward_sign``
        is not ``±1``. Both are fail-loud rather than clamped: a caller that
        does not know its own face's orientation cannot be handed a guess.

    Examples
    --------
    >>> face_name(0, -1), face_name(1, +1)
    ('xmin', 'ymax')
    """
    if not 0 <= axis < len(AXIS_NAMES):
        raise ValueError(
            f"face_name: axis index {axis} is outside the named-axis "
            f"inventory {AXIS_NAMES}."
        )
    try:
        suffix = _FACE_SIGN_SUFFIX[int(outward_sign)]
    except KeyError:
        raise ValueError(
            f"face_name: outward_sign must be +1 (the axis's 'max' face) or "
            f"-1 (its 'min' face); got {outward_sign!r}."
        ) from None
    return f"{AXIS_NAMES[axis]}{suffix}"


def face_normal(face: str) -> tuple[int, int]:
    r"""Parse a canonical face name into ``(axis index, outward-normal sign)``.

    The reverse direction of the face-name bijection: ``"xmin" -> (0, -1)``.
    Inverse of :func:`face_name`. The returned pair IS the face's outward
    normal :math:`\hat n_f = \mathrm{sign}\cdot\hat e_{\rm axis}`, stored
    sparsely — which is why
    :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`
    can derive the whole signed-projection table from face NAMES alone.

    Raises
    ------
    ValueError
        If ``face`` is not a canonical face name.

    Examples
    --------
    >>> face_normal("xmin"), face_normal("ymax")
    ((0, -1), (1, 1))
    """
    for suffix, sign in _FACE_SUFFIX_SIGN.items():
        if face.endswith(suffix):
            stem = face[: -len(suffix)]
            if stem in AXIS_NAMES:
                return AXIS_NAMES.index(stem), sign
    raise ValueError(
        f"Unknown face name {face!r}; valid faces are {sorted(FACE_NAMES)}"
    )


def face_opposite(face: str) -> str:
    r"""The face across the domain: same axis, opposite outward normal.

    ``"xmin" -> "xmax"``. An involution, and the only thing a *translation*
    deck transformation needs to know about where it was installed — which is
    why :class:`~orpheus.geometry.boundary.PairedDeck` carries the motion
    rather than a partner face (the partner is configuration, the motion is
    intrinsic).

    The two faces' outward normals are opposite (:math:`\hat n_f =
    -\hat n_{f'}`), so a direction OUTGOING at one is INCOMING at the other —
    the reason a periodic pair's crossing comes for free, with
    :math:`|\Omega\cdot\hat n|` preserved.

    .. note::

       This is the *geometric* opposite, not a gluing map. A non-opposite
       identification — a hex partner, a rotational quotient — is a different
       object and needs an explicit partner map (issue **#178**).

    Examples
    --------
    >>> face_opposite("xmin"), face_opposite(face_opposite("xmin"))
    ('xmax', 'xmin')
    """
    axis, sign = face_normal(face)
    return face_name(axis, -sign)


#: Every canonical face name, in axis-then-endpoint order. The inventory the
#: bijection is total on; consumers listing "the valid faces" read it rather
#: than re-deriving the product.
FACE_NAMES: tuple[str, ...] = tuple(
    face_name(axis, sign)
    for axis in range(len(AXIS_NAMES))
    for sign in (-1, +1)
)


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

    def pack(self, face_arrays: Mapping[K, NDArray]) -> NDArray:
        r"""Pack per-face ndarrays into a fresh flat backing buffer.

        ``face_arrays`` must cover EVERY face of this layout, each array
        matching its slot's shape. The packing loop lives HERE because
        every fact it reads — the face set, each slot's shape, the slice
        map, the total size — is layout knowledge (native place, CS4b
        S6.2); the field factories
        (:meth:`~orpheus.transport.fields._bases.FaceField.from_face_arrays`)
        are typed entries over this.

        Raises
        ------
        ValueError
            If the keys differ from this layout's faces, or any per-face
            array's shape mismatches its slot.
        """
        provided = set(face_arrays.keys())
        expected = set(self.faces.keys())
        if provided != expected:
            raise ValueError(
                f"FaceLayout.pack: face_arrays keys "
                f"{sorted(provided, key=repr)!r} do not match the "
                f"layout's faces {sorted(expected, key=repr)!r}"
            )
        flat = np.zeros(self.total_size)
        for name, slot in self.faces.items():
            arr = np.asarray(face_arrays[name])
            if arr.shape != slot.shape:
                raise ValueError(
                    f"FaceLayout.pack: face {name!r} array shape "
                    f"{arr.shape!r} does not match its slot shape "
                    f"{slot.shape!r}"
                )
            slot.slice_view(flat)[:] = arr
        return flat

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
