r"""The colouring symmetry of a layered geometry: the part of its mirror group that the colouring respects.

A 1-D layered geometry is an interval :math:`[b_0, b_n]` coloured by material runs and closed at its two ends
by face roles (a deck or a tag, :class:`~orpheus.geometry.deck.FaceRole`). The reflection
:math:`x \mapsto (b_0 + b_n) - x` is a symmetry of the *problem* when it maps the colouring and the closure
onto themselves; the group :math:`H` it generates (identity first, then the reflection when admitted) is the
stabiliser of the fundamental domain.

**When the reflection is in** :math:`H`. With the breakpoints of the *runs* :math:`b_0 < \dots < b_n`
(adjacent equal material ids merged), all three must hold:

1. every mirrored breakpoint :math:`m_i = (b_0 + b_n) - b_{n-i}` agrees with :math:`b_i` under the geometry's one
   coordinate-agreement convention, :func:`~orpheus.geometry.chart.coordinates_agree`, at the scale
   :math:`\max(|b_0|, |b_n|)` (64 ulp there; measured in ``measure_float_rule.md``: the exact comparison accepts
   2 of 13 symmetric slabs, a width-relative tolerance refuses an offset slab, 64 ulp at the end magnitude
   accepts 13 of 13 and refuses every asymmetric one, the genuine ceiling being 3.0 ulp);
2. the material runs read the same reversed, compared EXACTLY (integer ids);
3. the two end roles are the same tag or the same deck type.

**The runs are merged here, once.** The three producers of a raw partition (a geometry's intervals, a mesh's
regions, a material mesh's cells) pass what they hold and this function merges equal adjacent ids, so a group is
a function of the colouring and not of the cell partition: refining a mesh cannot change admission.

**What the value is.** For a slab whose two faces carry :class:`~orpheus.geometry.deck.TranslationDeck` the
orbifold is a circle, which has more reflections than :math:`x \mapsto a + b - x`; the value is therefore the
*stabiliser of the fundamental domain*, and stays that until a full deck group lands. A radial chart returns the
identity only (:math:`G_c` absorbs the rotation and the quotient is trivial). A multi-axis carrier has no group
machinery yet and answers :class:`UncomputedSymmetries`, the typed pessimistic answer, never an optimistic
:math:`\{e\}`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TypeAlias

import numpy as np

from orpheus.geometry.chart import coordinates_agree
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.deck import DeckTransformation, FaceRole
from orpheus.geometry.faces import face_inventory, face_name
from orpheus.geometry.transformation import RigidMotion
from orpheus.values.content import ContentIdentity, FrozenMapping

__all__ = [
    "ColouringSymmetry",
    "FaceIsometry",
    "Symmetries",
    "UncomputedSymmetries",
    "colouring_symmetry",
]


@dataclass(frozen=True, eq=False)
class FaceIsometry(ContentIdentity):
    r"""An affine isometry of the domain and the permutation it induces on the boundary faces.

    ``motion`` carries its position (:math:`x \mapsto a + b - x`); ``faces`` maps each boundary face to its image.
    """

    motion: RigidMotion
    faces: FrozenMapping[str, str]

    def __post_init__(self) -> None:
        if set(self.faces.values()) != set(self.faces):
            raise ValueError(f"FaceIsometry: the face map {dict(self.faces)} is not a permutation of its own faces")

    @property
    def is_identity(self) -> bool:
        """The motion is the identity and every face is fixed."""
        return self.motion.is_identity and all(k == v for k, v in self.faces.items())


@dataclass(frozen=True, eq=False)
class ColouringSymmetry(ContentIdentity):
    """The stabiliser of the fundamental domain: ``elements`` with the identity first (one or two for a slab)."""

    elements: tuple[FaceIsometry, ...]

    def __post_init__(self) -> None:
        if not self.elements or not self.elements[0].is_identity:
            raise ValueError("ColouringSymmetry: the elements are non-empty and the identity comes first")


@dataclass(frozen=True, eq=False)
class UncomputedSymmetries(ContentIdentity):
    """The typed pessimistic answer where no group machinery exists yet; ``carrier`` names the carrier, e.g. ``"Mesh2D"``."""

    carrier: str


Symmetries: TypeAlias = ColouringSymmetry | UncomputedSymmetries

_END_FACES = (face_name(0, -1), face_name(0, +1))


def _same_role(low: FaceRole, high: FaceRole) -> bool:
    """The same tag, or the same deck type; a deck never stands for a tag (a falsy tag is still a tag)."""
    if isinstance(low, DeckTransformation) or isinstance(high, DeckTransformation):
        return type(low) is type(high)
    return low == high


def colouring_symmetry(
    coord: CoordSystem,
    ends: tuple[float, float],
    breaks: Sequence[float],
    ids: Sequence[int],
    role_low: FaceRole,
    role_high: FaceRole,
) -> ColouringSymmetry:
    r"""The group of a 1-D colouring, from the RAW partition ``breaks`` (all :math:`n + 1` of them, ends included)
    and ``ids`` (one material id per interval); equal adjacent ids are merged into runs here.

    Parameters
    ----------
    coord:
        The coordinate system; only a Cartesian chart can have a mirror in :math:`H`.
    ends:
        ``(b_0, b_n)``, the interval's end points (the scale and the mirror's offset); they must equal the first
        and last break, and are refused otherwise (one quantity, one spelling: a disagreeing pair would silently
        change the scale and the mirror's offset).
    breaks, ids:
        The raw partition: a geometry's breakpoints and material ids, a mesh's region edges and ids, a material
        mesh's cell edges and ``mat_map`` values.
    role_low, role_high:
        The roles of the low and the high end face.

    Raises
    ------
    ValueError
        If the partition is not parsed: no interval, ``len(breaks) != len(ids) + 1``, breaks not strictly
        increasing, or ``ends`` not the first and last break. Checked before the chart is read.
    """
    # The raw partition is parsed first, whatever the chart: one more break than intervals, at least one interval,
    # strictly increasing, and the stated ends are its first and last break (one spelling of b_0 and b_n).
    if len(ids) < 1 or len(breaks) != len(ids) + 1:
        raise ValueError(f"colouring_symmetry: {len(breaks)} breaks for {len(ids)} intervals; expected at least one interval and one more break than ids")
    if not np.all(np.diff(np.asarray(breaks, dtype=float)) > 0.0):
        raise ValueError(f"colouring_symmetry: the breaks {list(breaks)} are not strictly increasing")
    low, high = float(ends[0]), float(ends[1])
    if (low, high) != (float(breaks[0]), float(breaks[-1])):
        raise ValueError(f"colouring_symmetry: ends {(low, high)} disagree with the first and last break {(float(breaks[0]), float(breaks[-1]))}")
    # The identity fixes the faces the carrier HAS: both ends on a slab and on a hollow radial body, the outer surface alone on a
    # solid one (its centre is an interior point, not a face), the one inventory rule of ``face_inventory``. A group acting on a
    # face the carrier does not have would be refused by ``Binding`` as another carrier's group.
    identity = FaceIsometry(RigidMotion.identity(1), FrozenMapping({face: face for face in face_inventory(coord, (low, high), 1)}))
    if coord is not CoordSystem.CARTESIAN:
        return ColouringSymmetry((identity,))

    run_ids = [ids[0]]
    run_breaks = [float(breaks[0])]
    for i in range(1, len(ids)):
        if ids[i] != ids[i - 1]:
            run_breaks.append(float(breaks[i]))
            run_ids.append(ids[i])
    run_breaks.append(float(breaks[-1]))

    total = low + high
    scale = max(abs(low), abs(high))
    mirrored = [total - b for b in reversed(run_breaks)]
    symmetric = (
        bool(np.all(coordinates_agree(np.asarray(mirrored), np.asarray(run_breaks), scale=scale)))
        and all(a == b for a, b in zip(run_ids, reversed(run_ids)))
        and _same_role(role_low, role_high)
    )
    if not symmetric:
        return ColouringSymmetry((identity,))
    mirror = FaceIsometry(
        RigidMotion.reflection(normal=np.array([1.0]), offset=0.5 * total),
        FrozenMapping({_END_FACES[0]: _END_FACES[1], _END_FACES[1]: _END_FACES[0]}),
    )
    return ColouringSymmetry((identity, mirror))
