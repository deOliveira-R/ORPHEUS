r"""The deck transformations: how a domain continues across one of its faces.

A transport domain :math:`D` closes its boundary in exactly two ways (``wave2_boundary_carve.md`` section 1.1).
Either the domain *ends* at the face and an exterior model supplies the inflow (a boundary condition, a response
kernel plus a source, which is the transport layer's), or the domain *continues* across the face: the inflow
trace at :math:`f` is the outflow trace at a partner face :math:`f'` carried by an isometry :math:`g` of the
covering space,

.. math::

    \gamma_-\psi|_f = G_g\,\gamma_+\psi|_{f'},

with :math:`G_g` the Koopman operator of :math:`g`. :math:`G_g` is invertible, measure preserving and has no
parameter beyond :math:`g`: it has no amplitude, no source, nothing to scale. This module is the geometry's half:
the isometry classes, the two decks that are minted, the pairing a realizer reads of a deck, and the parse of a
face role. It names no response and imports no transport vocabulary.

**The classes of isometry (Chasles).** A non-identity isometry of :math:`E^n` (:math:`n \le 3`) is in exactly one
of six conjugacy classes: translation, reflection, rotation, glide reflection, screw motion, rotoreflection.
Dimension 1 has two (translation, reflection), dimension 2 four (adding rotation and glide reflection),
dimension 3 six. :class:`IsometryClass` enumerates them and :data:`ISOMETRY_CLASSES_BY_DIMENSION` records which
exist where.

**Why only two are minted.** The admissible decks of a chart are :math:`N(G_c)/G_c`, the normaliser of the
chart's own symmetry group modulo that group; on a Cartesian axis that is the mirror and the translation, and on
a 1-D radial chart it is trivial (:func:`admitted_decks`: a reflecting curved wall is a response, never a deck).
A rotation's partner face is not derivable from the face (a quarter turn about a corner pairs ``xmin`` with
``ymin`` or ``ymax``), a glide reflection and a screw carry a translation part and a rotoreflection an order, so
their fields are designed with the 2-D and 3-D geometry values. Minting them now would mint four types with
guessed fields and four unreachable arms in every method. The classification stays complete because the
deferred classes are rows of :data:`DEFERRED_DECKS`, keyed by the Enum, each naming the type it will have and the
machinery that mints it.

**A deck has no field.** The face fixes the plane (mirror), the partner and the direction of the translation, so
:class:`MirrorDeck` and :class:`TranslationDeck` store nothing and compare equal to themselves by type. A wrong
axis declaration (the class of defect where a mirror carried an axis different from its face's) is unspellable.
The position of a mirror plane is not stored: a 1-D mirror acts on directions only, because the offset is
invisible to every discrete quantity.

**Dispatch.** :data:`DeckTransformation` is a closed sum type and every method that reads one does so by an
exhaustive ``match`` ending in ``case _: assert_never(...)``; there is no registry. The day a third deck is
minted, every such match fails the type checker until it decides.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import ClassVar, Mapping, TypeAlias, assert_never

import numpy as np

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.faces import AXIS_INDEX, AXIS_NAMES, face_normal, face_opposite
from orpheus.geometry.transformation import RigidMotion
from orpheus.values.content import ContentIdentity

__all__ = [
    "DEFERRED_DECKS",
    "ISOMETRY_CLASSES_BY_DIMENSION",
    "DeckError",
    "DeckTransformation",
    "DeferredDeck",
    "FacePairing",
    "FaceRole",
    "IsometryClass",
    "MirrorDeck",
    "TranslationDeck",
    "admitted_decks",
    "check_pairings_close",
    "mirror_motion",
    "pairing_of",
    "parse_face_role",
    "translation_motion",
]


class DeckError(ValueError):
    """The geometry's refusal of an ill-posed deck declaration (not a boundary-condition error)."""


class IsometryClass(Enum):
    """The six conjugacy classes of non-identity isometries of :math:`E^n`, :math:`n \\le 3` (Chasles)."""

    TRANSLATION = "translation"
    REFLECTION = "reflection"
    ROTATION = "rotation"
    GLIDE_REFLECTION = "glide reflection"
    SCREW_MOTION = "screw motion"
    ROTOREFLECTION = "rotoreflection"


#: The classes that exist in each dimension; the classes of dimension :math:`n` are a proper subset of those of
#: :math:`n + 1`.
ISOMETRY_CLASSES_BY_DIMENSION: Mapping[int, frozenset[IsometryClass]] = MappingProxyType(
    {
        1: frozenset({IsometryClass.TRANSLATION, IsometryClass.REFLECTION}),
        2: frozenset(
            {IsometryClass.TRANSLATION, IsometryClass.REFLECTION, IsometryClass.ROTATION, IsometryClass.GLIDE_REFLECTION}
        ),
        3: frozenset(IsometryClass),
    }
)


@dataclass(frozen=True, eq=False)
class MirrorDeck(ContentIdentity):
    r"""The face is glued to itself by the reflection in its own plane: the domain continues as its mirror image.

    No field: the face fixes the plane. Two mirror decks are equal whichever face they sit on, and generate the
    same group :math:`\mathbb{Z}_2`.
    """

    isometry_class: ClassVar[IsometryClass] = IsometryClass.REFLECTION


@dataclass(frozen=True, eq=False)
class TranslationDeck(ContentIdentity):
    r"""The face is glued to the opposite face of the same axis by a translation: the domain is periodic.

    No field: the face fixes the partner (:func:`~orpheus.geometry.faces.face_opposite`) and the direction. The
    length of the translation is not content (a deck does not know the domain's extent).
    """

    isometry_class: ClassVar[IsometryClass] = IsometryClass.TRANSLATION


#: The minted decks, a closed sum.
DeckTransformation: TypeAlias = MirrorDeck | TranslationDeck

#: What a boundary face carries: a deck (the domain continues) or an ``int`` tag (the domain ends; the tag names
#: the exterior model the problem binds to it). A ``bool`` is not a tag.
FaceRole: TypeAlias = DeckTransformation | int


@dataclass(frozen=True)
class DeferredDeck:
    """An isometry class named and not minted: the type it will have and the machinery that mints it."""

    isometry_class: IsometryClass
    name: str
    trigger: str


#: The four classes the minted decks do not cover. Keyed by the Enum, so the classification is complete:
#: ``minted | deferred == set(IsometryClass)``.
DEFERRED_DECKS: Mapping[IsometryClass, DeferredDeck] = MappingProxyType(
    {
        IsometryClass.ROTATION: DeferredDeck(
            IsometryClass.ROTATION, "RotationDeck", "the 2-D geometry value (a rotation's partner face is not derivable from the face)"
        ),
        IsometryClass.GLIDE_REFLECTION: DeferredDeck(
            IsometryClass.GLIDE_REFLECTION, "GlideDeck", "the 2-D geometry value (a glide carries a translation part)"
        ),
        IsometryClass.SCREW_MOTION: DeferredDeck(
            IsometryClass.SCREW_MOTION, "ScrewDeck", "a 3-D geometry value (a screw carries an axis and a translation part)"
        ),
        IsometryClass.ROTOREFLECTION: DeferredDeck(
            IsometryClass.ROTOREFLECTION, "RotoreflectionDeck", "a 3-D geometry value (a rotoreflection carries an order)"
        ),
    }
)


def _axis_index(axis: str, dimension: int, what: str) -> int:
    try:
        index = AXIS_INDEX[axis]
    except KeyError:
        raise ValueError(f"{what}: axis must be one of {AXIS_NAMES}, got {axis!r}") from None
    if not 0 <= index < dimension:
        raise ValueError(f"{what}: axis {axis!r} is index {index}, out of range for dimension {dimension}")
    return index


def mirror_motion(axis: str, dimension: int) -> RigidMotion:
    r"""The reflection whose plane is normal to the named ``axis``: the ONE axis-letter to mirror spelling.

    The axis letter is resolved through :data:`~orpheus.geometry.faces.AXIS_INDEX` and refused here, so a
    mis-named axis cannot be constructed and fail later holding a quadrature. The plane passes through the
    origin: its position is not content (see the module docstring).

    Raises
    ------
    ValueError
        If ``axis`` is not a letter of :data:`~orpheus.geometry.faces.AXIS_NAMES` or its index is not below
        ``dimension``.
    """
    index = _axis_index(axis, dimension, "mirror_motion")
    normal = np.zeros(dimension)
    normal[index] = 1.0
    return RigidMotion.reflection(normal=normal)


def translation_motion(face: str, dimension: int) -> RigidMotion:
    r"""The unit translation along the outward normal of ``face``: :math:`g_f = \operatorname{sign}\,\hat e_{\rm axis}`.

    SIGNED BY THE FACE. The motion :math:`g_f` carries the partner face onto ``face``, so ``g_xmin`` is
    :math:`-\hat e_x` and ``g_xmax`` is :math:`+\hat e_x`, and the groupoid law :math:`g_{f'} = g_f^{-1}` closes.
    The vector's length is not content: readers use ``linear`` and ``is_translation`` only.

    Raises
    ------
    ValueError
        If ``face`` is not a canonical face name or its axis index is not below ``dimension``.
    """
    index, sign = face_normal(face)
    _axis_index(AXIS_NAMES[index], dimension, "translation_motion")
    vector = np.zeros(dimension)
    vector[index] = float(sign)
    return RigidMotion.translation_by(vector)


@dataclass(frozen=True)
class FacePairing:
    r"""What a realizer reads of a deck: the face, its partner and the rigid motion carrying the partner onto it.

    Derived from a deck and a face by :func:`pairing_of`, never stored on the deck. A face paired with itself
    (``partner == face``) carries an involution (a mirror).

    Raises
    ------
    DeckError
        If the face is its own partner and the motion's linear part is not an involution.
    """

    face: str
    partner: str
    motion: RigidMotion

    def __post_init__(self) -> None:
        if self.is_self_paired:
            linear = self.motion.linear
            if not np.array_equal(linear @ linear, np.eye(linear.shape[0])):
                raise DeckError(f"face {self.face!r} is paired with itself by a motion whose linear part is not an involution")

    @property
    def is_self_paired(self) -> bool:
        """The face is its own partner (a mirror)."""
        return self.partner == self.face

    @property
    def permutes_ordinates(self) -> bool:
        """The motion moves directions (it is not a translation): the single source is ``RigidMotion.is_translation``."""
        return not self.motion.is_translation

    @property
    def inverse(self) -> "FacePairing":
        r"""The pairing seen from the partner: :math:`(f', f, g^{-1})`, the adjoint's pairing (law L-A)."""
        return FacePairing(self.partner, self.face, self.motion.inverse())


def pairing_of(deck: DeckTransformation, face: str, *, dimension: int) -> FacePairing:
    """The pairing of ``face`` under ``deck`` in a ``dimension``-dimensional space (one exhaustive match)."""
    match deck:
        case MirrorDeck():
            axis, _ = face_normal(face)
            return FacePairing(face, face, mirror_motion(AXIS_NAMES[axis], dimension))
        case TranslationDeck():
            return FacePairing(face, face_opposite(face), translation_motion(face, dimension))
        case _:
            assert_never(deck)


def admitted_decks(coord: CoordSystem, axis: int) -> frozenset[type[DeckTransformation]]:
    r"""The deck types a face of ``axis`` admits on ``coord``: :math:`N(G_c)/G_c`.

    A Cartesian axis, and any non-radial axis of a curvilinear chart, admits the mirror and the translation. The
    radial first axis of a cylinder or a sphere admits none: :math:`G_c` absorbs the rotational symmetry and the
    quotient is trivial, so a reflecting curved wall is a response (a tagged face with a specular wall
    condition), never a deck.
    """
    if coord is not CoordSystem.CARTESIAN and axis == 0:
        return frozenset()
    return frozenset({MirrorDeck, TranslationDeck})


def parse_face_role(role: object, where: str) -> FaceRole:
    """A deck, or an ``int`` tag (numpy integers coerced to ``int``); anything else is a ``TypeError`` naming ``where``.

    A ``bool`` is not a tag (the twin of the material-id rule, ``orpheus/data/materials.py``). No sign restriction.
    """
    if isinstance(role, MirrorDeck | TranslationDeck):
        return role
    if isinstance(role, int | np.integer) and not isinstance(role, bool):
        return int(role)
    raise TypeError(f"{where}: a face role is a MirrorDeck, a TranslationDeck or an int tag, got {role!r} of type {type(role).__name__}")


def check_pairings_close(roles: Mapping[str, FaceRole], where: str) -> None:
    """A translation on a face requires a translation on the opposite face (an involution); else ``DeckError``."""
    for face, role in roles.items():
        if isinstance(role, TranslationDeck):
            partner = face_opposite(face)
            if not isinstance(roles.get(partner), TranslationDeck):
                raise DeckError(f"{where}: face {face!r} carries a TranslationDeck but its partner {partner!r} carries {roles.get(partner)!r}")
