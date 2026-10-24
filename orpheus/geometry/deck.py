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

**Dispatch.** :data:`DeckTransformation` is a closed sum type and there is no registry. Two mechanisms read it,
and they govern different things. BEHAVIOUR (what a deck does: its pairing) is an exhaustive ``match`` ending in
``case _: assert_never(...)`` (:func:`pairing_of`), so the day a third deck is minted the type checker fails that
match until it decides. MEMBERSHIP (is this value a deck, which deck types does a chart admit) reads the alias
itself, by ``isinstance(x, DeckTransformation)`` and ``get_args(DeckTransformation)``, so a minted deck is a member
the moment the alias names it and no site re-spells the union by hand.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from collections.abc import Iterable
from typing import ClassVar, Mapping, TypeAlias, assert_never, get_args

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
    "parse_face_role_on",
    "parse_face_roles",
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
    r"""The reflection whose plane is normal to the named ``axis``: the geometry's axis-letter to mirror spelling.

    The one axis-letter to mirror spelling: ``orpheus/geometry/symmetry.py`` (``_reflections``) delegates here, keeping its
    message's fragment "axis must be one of" (``wave2_boundary_carve.md`` section 5.1; wave-2 slice S1; the second
    delegate, ``_mirror_motion`` of the retired ``orpheus/geometry/boundary/_factors.py``, retired in WF2.4).

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

    Derived from a deck and a face by :func:`pairing_of`, never stored on the deck. The class states the whole law,
    so a hand-built value cannot spell a defect :func:`pairing_of` cannot:

    * both faces are canonical face names (:func:`~orpheus.geometry.faces.face_normal`);
    * a face paired with itself carries EXACTLY the reflection in its own plane,
      ``mirror_motion(axis of the face, dimension)`` (neither the identity, nor a translation, nor a mirror of
      another axis);
    * a pair of distinct faces whose motion is a translation (linear part the identity) is the periodic gluing:
      the partner is ``face_opposite(face)`` and the translation is non-zero, along the face's axis, with the
      face's outward sign;
    * a pair of distinct faces whose motion is not a translation (a rotation, a glide, a screw or a
      rotoreflection) is refused: it belongs to a class of :data:`DEFERRED_DECKS`, and a pairing exists only for
      a minted deck.

    SCOPE-BOUNDARY[guard] machinery: the deferred decks of :data:`DEFERRED_DECKS`, each stating its own pairing law here.
    ruling: the orchestrator, wave 2 WF2.2 integration, 2026-10-23 (a pairing exists only for a minted deck; a quarter turn glues adjacent faces, so no law of today's admits it).
    revisit: when the first deferred deck is minted, its row adds its arm to this law.

    Raises
    ------
    DeckError
        If a face name is not canonical, or the law above fails.
    """

    face: str
    partner: str
    motion: RigidMotion

    def __post_init__(self) -> None:
        try:
            axis, sign = face_normal(self.face)
            face_normal(self.partner)
        except ValueError as error:
            raise DeckError(f"a pairing names a face that is not canonical: {error}") from None
        dimension = self.motion.linear.shape[0]
        if axis >= dimension:
            raise DeckError(f"face {self.face!r} is on axis {axis}, out of range for a motion of dimension {dimension}")
        if self.is_self_paired:
            expected = mirror_motion(AXIS_NAMES[axis], dimension)
            if not (
                np.array_equal(self.motion.linear, expected.linear)
                and np.array_equal(self.motion.translation, expected.translation)
            ):
                raise DeckError(f"face {self.face!r} is paired with itself by a motion that is not the reflection in its own plane")
        elif self.motion.is_translation:
            direction = np.zeros(dimension)
            direction[axis] = float(sign)
            along = self.motion.translation
            if self.partner != face_opposite(self.face) or not (
                np.any(along != 0.0) and np.array_equal(np.sign(along), direction)
            ):
                raise DeckError(
                    f"face {self.face!r} is glued to {self.partner!r} by a translation {along.tolist()}: a translation pairs a face "
                    f"with its opposite, along its axis, with its outward sign"
                )
        else:
            raise DeckError(
                f"face {self.face!r} is glued to {self.partner!r} by a motion that is neither its own mirror nor a translation: "
                f"that is a deferred deck (DEFERRED_DECKS), and a pairing exists only for a minted deck"
            )

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
    return frozenset(get_args(DeckTransformation))


def parse_face_role(role: object, where: str) -> FaceRole:
    """A deck, or an ``int`` tag (numpy integers coerced to ``int``); anything else is a ``TypeError`` naming ``where``.

    A ``bool`` is not a tag (the twin of the material-id rule, ``orpheus/data/materials.py``). No sign restriction.
    """
    if isinstance(role, DeckTransformation):
        return role
    if isinstance(role, int | np.integer) and not isinstance(role, bool):
        return int(role)
    raise TypeError(f"{where}: a face role is a MirrorDeck, a TranslationDeck or an int tag, got {role!r} of type {type(role).__name__}")


def check_pairings_close(roles: Mapping[str, FaceRole], where: str) -> None:
    """A deck gluing a face to a DIFFERENT face requires the same deck type there (an involution); else ``DeckError``.

    The partner is read from :func:`pairing_of`, never re-derived (today only the translation glues two faces; a
    mirror is self-paired and has nothing to close). ``pairing_of`` is asked in the smallest dimension that holds
    the face's axis: the partner does not depend on the dimension.
    """
    for face, role in roles.items():
        if isinstance(role, DeckTransformation):
            partner = pairing_of(role, face, dimension=face_normal(face)[0] + 1).partner
            if partner != face and type(roles.get(partner)) is not type(role):
                raise DeckError(f"{where}: face {face!r} carries a {type(role).__name__} but its partner {partner!r} carries {roles.get(partner)!r}")


def tags_of(roles: Iterable[FaceRole]) -> frozenset[int]:
    """The tags among ``roles``, the decks excluded: the ONE predicate behind every ``tags`` property."""
    return frozenset(role for role in roles if not isinstance(role, DeckTransformation))


def parse_face_role_on(role: object, coord: CoordSystem, axis: int, where: str) -> FaceRole:
    r"""One face role, parsed and admitted on the ``axis`` of ``coord``: :func:`parse_face_role`, then the chart's deck rule.

    A deck on a face whose axis admits none (:func:`admitted_decks`: the radial first axis of a cylinder or a sphere)
    is a :class:`DeckError` naming the quotient :math:`W = N(G_c)/G_c`, which is trivial there. The one place the
    refusal is spelled: the geometry, both meshes and both axis classes read their roles through it.
    """
    parsed = parse_face_role(role, where)
    if isinstance(parsed, DeckTransformation) and type(parsed) not in admitted_decks(coord, axis):
        raise DeckError(
            f"{where}: a {type(parsed).__name__} on axis {axis} of a {coord.name.lower()} chart: "
            f"W = N(G_c)/G_c is trivial on the radial axis of a curvilinear chart: a reflecting curved wall is a response; "
            f"bind a specular exterior condition to a tag"
        )
    return parsed


def parse_face_roles(roles: Mapping[str, object], coord: CoordSystem, where: str) -> dict[str, FaceRole]:
    """The roles of named faces, each through :func:`parse_face_role_on` at its own axis, then the pairings closed once.

    ``roles`` maps a canonical face name to its declared role; the order of the result is the order given. The count
    of faces against the chart's inventory is the caller's (the geometry's boundary points, the mesh's
    ``face_inventory``): this function reads only what each face and its partner carry.
    """
    parsed = {face: parse_face_role_on(role, coord, face_normal(face)[0], f"{where}[{face!r}]") for face, role in roles.items()}
    check_pairings_close(parsed, where)
    return parsed
