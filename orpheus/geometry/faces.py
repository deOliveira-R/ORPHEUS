r"""The boundary-face naming of an axis-aligned domain: the axis letters and the face bijection.

A boundary face of an axis-aligned domain is completely described by its
OUTWARD NORMAL, and for an axis-aligned face that normal is sparse,
:math:`\hat n_f = \mathrm{sign}\cdot\hat e_{\rm axis}`. This module holds the
two encodings of that object and the bijection between them: the pair
``(axis index, outward sign)`` and the face NAME ``"{axis}{min|max}"``.

It is geometry (layer DAG WP-B, 2026-10-22, ruling 1 of
``.claude/plans/layer_dag.md``): a face, its normal and the name of a
coordinate axis are statements about the domain, and need no measure. Until
2026-10-22 these names lived in ``orpheus/numerics/face_layout.py`` beside the
flat-buffer descriptor and the streaming-normal measure, which stay in numerics
(:mod:`orpheus.numerics.face_layout`) and import the names from here.

:data:`AXIS_NAMES` is the ONE spelling of :math:`\mathbb{R}^3`'s axis letters
in the tree. The point groups (:mod:`orpheus.geometry.symmetry`) read it to name
the axis of :math:`SO(2)_a` and :math:`O(2)_a` and the normal of
:math:`\sigma_a`; the face names read it for their stem. Until WP-B the point
groups held a second table (``AXIS_LETTER``, with its inverse ``AXIS_INDEX``),
equal to this one by transcription; the merge left one table and its derived
inverse.

This module imports nothing from orpheus.
"""

from __future__ import annotations

__all__ = [
    "AXIS_INDEX",
    "AXIS_NAMES",
    "FACE_NAMES",
    "face_name",
    "face_normal",
    "face_opposite",
]


#: The letters of :math:`\mathbb{R}^3`'s coordinate axes, positional-by-axis:
#: the single table of the axis↔letter crosswalk. Two conventions key on it.
#: The face names ``"{axis}min"`` / ``"{axis}max"`` that
#: :class:`~orpheus.numerics.face_layout.FaceLayout`, the trace space,
#: :attr:`SNProblem.bc <orpheus.sn.problem.SNProblem.bc>` and the sweep
#: schedule key on render through it, and no consumer hand-lists
#: ``("x", ...), ("y", ...)`` pairs (every face-name derivation,
#: :attr:`FaceLabel.face_name <orpheus.mesh.axis.FaceLabel.face_name>`, the
#: walk's in/outflow faces, the schedule's outgoing faces, the trace's
#: outward-normal table, renders through this tuple). The point groups of
#: :mod:`orpheus.geometry.symmetry` name an axis by its letter, every
#: realization built in this frame (principal axis along z, vertex on x); the
#: mirror family names its plane by the NORMAL, so ``AXIS_INDEX["x"]`` selects
#: the normal of :math:`\sigma_x` as well as the axis of :math:`SO(2)_x` and
#: :math:`O(2)_x`.
AXIS_NAMES: tuple[str, str, str] = ("x", "y", "z")

#: The inverse, ``AXIS_INDEX[AXIS_NAMES[i]] == i``, DERIVED from the letters so
#: the two directions cannot disagree.
AXIS_INDEX: dict[str, int] = {a: i for i, a in enumerate(AXIS_NAMES)}


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
# It lives beside :data:`AXIS_NAMES` because the face name is the axis letter
# plus the sign's suffix (the trace space and :mod:`orpheus.mesh.axis` both
# import it). Before campaign phase **B3.4c** the convention had no single
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
